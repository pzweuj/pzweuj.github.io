# OPDS/Mihon 读取 EPUB 章节返回 404

## 版本

Kavita v0.9.0.2（`jvmilazz0/kavita:latest`）

## 现象

通过 OPDS（Mihon 扩展）打开 EPUB 格式的漫画/书籍时，页面图片请求返回 404。

Mihon 发出的请求：
```
GET /api/Reader/image?chapterId=3&page=0&extractPdf=true&apiKey=...
```

Kavita 日志：
```
[Information] Serilog.AspNetCore.RequestLoggingMiddleware
HTTP GET "/api/Reader/image?chapterId=3&page=0&extractPdf=true&apiKey=******REDACTED******"
responded 404 in 3.8976 ms
```

同一章节通过 Kavita Web 阅读器可以正常打开，`/api/Reader/pdf` 可以正常下载 EPUB 文件。

## 根因

Kavita **不支持**服务端 EPUB → 图片的逐页提取。

`ReadingItemService.Extract()` 对 EPUB 格式的处理是**空操作**：

**文件：** `Kavita.Services/Reading/ReadingItemService.cs`

```csharp
// v0.9.0.2 tag
public void Extract(string fileFilePath, string targetDirectory, MangaFormat format, int imageCount = 0)
{
    switch (format)
    {
        case MangaFormat.Archive:
            _archiveService.ExtractArchive(fileFilePath, targetDirectory);
            break;

        case MangaFormat.Image:
            _imageService.ExtractImages(fileFilePath, targetDirectory, imageCount);
            break;

        case MangaFormat.Pdf:
            _bookService.ExtractPdfImages(fileFilePath, targetDirectory);
            break;

        case MangaFormat.Epub:
            break;  // <-- 空操作，未实现 EPUB 页面提取

        case MangaFormat.Unknown:
        default:
            break;
    }
}
```

## 调用链

```
Mihon Extension (Kavita.kt)
  │
  │  构造 URL: $apiUrl/Reader/image?chapterId=3&page=0&extractPdf=true&apiKey=$apiKey
  │
  ▼
ReaderController.GetImage()
  │  Kavita.Server/Controllers/ReaderController.cs  [HttpGet("image")]
  │
  ├─► [ChapterAccess]           → Admin 用户直接绕过
  ├─► [SkipDeviceTracking]
  │
  └─► cacheService.Ensure(chapterId, extractPdf=true)
        │  Kavita.Services/CacheService.cs
        │
        └─► ExtractChapterFiles(extractPath, files, extractPdfToImages=true)
              │
              └─► readingItemService.Extract(filePath, dir, MangaFormat.Epub)
                    │  Kavita.Services/Reading/ReadingItemService.cs
                    │
                    └─► case MangaFormat.Epub: break;  // 什么都不做
                          │
                          ▼
                    缓存目录只存在 EPUB 原文件，无页面图片
                          │
                          ▼
              GetCachedPagePath(chapter.Id, page)
                │  Kavita.Services/CacheService.cs
                │
                │  扫描缓存目录中的图片文件 → 找不到 → 返回 string.Empty
                │
                ▼
              CachedFile("")
                │  Kavita.Server/Controllers/BaseApiController.cs
                │
                │  if (string.IsNullOrEmpty(path) || !File.Exists(path))
                │      return NotFound();  // 404
                │
                ▼
              HTTP 404
```

## 为什么 Web 阅读器正常

Kavita Web 阅读器处理 EPUB 时，是将 EPUB 内的 HTML/XHTML 内容直接发送到浏览器，由浏览器端 JavaScript（epub.js 或类似库）渲染。不需要服务端提取图片。

`BookService`（`Kavita.Services/BookService.cs`）的作用是：
- 解析 EPUB 的 OPF 文件获取阅读顺序和页数
- 重写 HTML 中的资源引用（CSS、图片等）使其指向 API 端点
- 建立内容文件到页码的映射

这些操作不产生图片文件。Web 阅读器通过 `/api/Reader/pdf` 下载原始 EPUB，在浏览器中逐页渲染 HTML。

## 各格式 `/api/Reader/image` 支持情况

| MangaFormat | 格式 | 提取方式 | 支持 |
|-------------|------|---------|------|
| `Archive` | CBZ/CBR/ZIP/7Z | 解压提取内部图片文件 | ✅ |
| `Pdf` | PDF | 通过 `Docnet.Core` 渲染页面 | ✅ |
| `Image` | 单张图片 | 直接复制文件 | ✅ |
| `Epub` | EPUB | `break;`（空操作） | ❌ |
| `Unknown` | 未知 | `break;`（空操作） | ❌ |

## Mihon 扩展侧

**文件：** `tach-extension/src/all/kavita/src/eu/kanade/tachiyomi/extension/all/kavita/Kavita.kt`

扩展在 `getPageList()` 中为所有格式统一构造 `/api/Reader/image` URL：

```kotlin
// 不区分格式，统一使用 Reader/image
(0 until chapterDetails.pages).map { i ->
    Page(
        index = i,
        imageUrl = "$apiUrl/Reader/image?chapterId=$chapterId&page=$i&extractPdf=true&apiKey=$apiKey",
    )
}
```

扩展没有针对 EPUB 格式做特殊处理（如通过 OPDS 端点或下载原文件）。这导致 EPUB 章节的每页请求都是 404。

## 解决方案

1. **将 EPUB 转换为 CBZ**（推荐）
   ```bash
   ebook-convert "input.epub" "output.cbz"
   ```
   转换后重新导入 Kavita。CBZ 格式 `/api/Reader/image` 完全支持。

2. **EPUB 用 Kavita Web 阅读器**，CBZ/PDF 用 Mihon。

3. **向 Kavita 提 Issue**，建议实现 `ReadingItemService.Extract` 的 EPUB 分支，可通过集成 calibre CLI 或 headless browser 逐页渲染 EPUB 为图片。

4. **向 tach-extension 提 Issue**，建议对 EPUB 格式降级为下载整个 EPUB 文件（`/api/Reader/pdf`），由 Mihon 本地渲染。
