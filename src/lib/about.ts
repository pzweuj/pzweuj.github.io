import fs from 'fs'
import path from 'path'
import matter from 'gray-matter'
import { createMarkdownProcessor, renderMarkdown, PROCESSOR_VERSION } from './createMarkdownProcessor'
import { loadHtmlCache, saveHtmlCache, renderMarkdownWithCache } from './cache'

// 关于页管道：基础渲染
const processor = createMarkdownProcessor()

let aboutCache: { title: string; content: string } | null = null

export async function getAboutContent() {
  if (aboutCache) return aboutCache

  const filePath = path.join(process.cwd(), 'content/about.md')
  const htmlCache = loadHtmlCache()
  const cacheKey = 'content/about.md'
  const content = fs.readFileSync(filePath, 'utf-8')
  const { data, content: markdown } = matter(content)

  const result = {
    title: data.title || '关于',
    content: await renderMarkdownWithCache(
      markdown, cacheKey, 'about', PROCESSOR_VERSION, htmlCache,
      (c) => renderMarkdown(processor, c, cacheKey),
    )
  }
  saveHtmlCache(htmlCache)
  aboutCache = result
  return result
}
