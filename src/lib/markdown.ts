import fs from 'fs'
import path from 'path'
import matter from 'gray-matter'
import { createMarkdownProcessor, renderMarkdown, PROCESSOR_VERSION } from './createMarkdownProcessor'
import { loadHtmlCache, saveHtmlCache, renderMarkdownWithCache } from './cache'

// 博客管道：启用表格标签、代码复制按钮、图片懒加载
const processor = createMarkdownProcessor({
  enableTableLabel: true,
  enableCodeCopy: true,
  enableLazyImage: true,
})

export interface BlogPost {
  slug: string
  title: string
  date: string
  tags: string[]
  excerpt: string
  content: string
}

// 从文件名中提取日期和slug
function parseFileName(fileName: string) {
  const match = fileName.match(/^(\d{4}-\d{2}-\d{2})-(.+)\.md$/)
  if (!match) throw new Error(`Invalid file name format: ${fileName}`)
  return {
    date: match[1],
    slug: match[2]
  }
}

// 生成摘要
function generateExcerpt(content: string, maxLength: number = 200) {
  const plainText = content
    .replace(/<[^>]+>/g, '')
    .replace(/[#*`]/g, '')
    .trim()

  const firstParagraph = plainText
    .split('\n')
    .find(line => line.trim().length > 0)
    ?.trim() || ''

  return firstParagraph.length > maxLength
    ? `${firstParagraph.slice(0, maxLength)}...`
    : firstParagraph
}

function parseTags(tags: unknown): string[] {
  if (typeof tags === 'string') {
    return tags
      .split(/[,\s]+/)
      .map(tag => tag.trim())
      .filter(Boolean)
  }

  if (Array.isArray(tags)) {
    return tags
      .map(tag => String(tag).trim())
      .filter(Boolean)
  }

  return []
}

// 内存缓存：避免同一次构建中重复调用
const postsMemoCache = new Map<string, BlogPost[]>()

export async function getAllPosts(): Promise<BlogPost[]> {
  const key = 'content/posts'
  if (postsMemoCache.has(key)) return postsMemoCache.get(key)!
  const posts = await getAllPostsFromDirectory(key)
  postsMemoCache.set(key, posts)
  return posts
}

// 通用的从指定目录获取所有文章的函数
export async function getAllPostsFromDirectory(directory: string): Promise<BlogPost[]> {
  const postsDirectory = path.join(process.cwd(), directory)

  if (!fs.existsSync(postsDirectory)) {
    return []
  }

  const htmlCache = loadHtmlCache()
  const fileNames = fs.readdirSync(postsDirectory)

  const posts = await Promise.all(fileNames
    .filter(fileName => fileName.endsWith('.md'))
    .map(async fileName => {
      const { date, slug } = parseFileName(fileName)
      const fullPath = path.join(postsDirectory, fileName)
      const cacheKey = path.relative(process.cwd(), fullPath).replace(/\\/g, '/')
      const fileContents = fs.readFileSync(fullPath, 'utf8')
      const { data, content } = matter(fileContents)

      return {
        slug,
        title: data.title || slug,
        date,
        tags: parseTags(data.tags),
        excerpt: generateExcerpt(content),
        content: await renderMarkdownWithCache(
          content, cacheKey, 'blog', PROCESSOR_VERSION, htmlCache,
          (c) => renderMarkdown(processor, c, cacheKey),
        )
      }
    }))

  saveHtmlCache(htmlCache)
  return posts.sort((a, b) => b.date.localeCompare(a.date))
}

// 获取 Schema 进度文章
export async function getAllSchemaProgress(): Promise<BlogPost[]> {
  const key = 'content/schema-progress'
  if (postsMemoCache.has(key)) return postsMemoCache.get(key)!
  const posts = await getAllPostsFromDirectory(key)
  postsMemoCache.set(key, posts)
  return posts
}

// 分页相关接口定义
export interface PaginationInfo {
  currentPage: number
  totalPages: number
  hasNextPage: boolean
  hasPrevPage: boolean
}

export interface PaginatedPosts {
  posts: BlogPost[]
  pagination: PaginationInfo
}

export async function getPaginatedPosts(page: number = 1, limit: number = 5): Promise<PaginatedPosts> {
  const posts = await getAllPosts()
  const startIndex = (page - 1) * limit
  const endIndex = startIndex + limit
  const totalPages = Math.ceil(posts.length / limit)

  return {
    posts: posts.slice(startIndex, endIndex),
    pagination: {
      currentPage: page,
      totalPages,
      hasNextPage: endIndex < posts.length,
      hasPrevPage: page > 1
    }
  }
}
