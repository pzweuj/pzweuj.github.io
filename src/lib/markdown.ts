import fs from 'fs'
import path from 'path'
import matter from 'gray-matter'
import { unified } from 'unified'
import remarkParse from 'remark-parse'
import remarkMath from 'remark-math'
import remarkGfm from 'remark-gfm'
import remarkRehype from 'remark-rehype'
import rehypeKatex from 'rehype-katex'
import rehypeStringify from 'rehype-stringify'
import rehypePrism from 'rehype-prism-plus'
import rehypeImgSize from 'rehype-img-size'
import { remarkQQMusic } from './plugins/remarkQQMusic'
import { rehypeTableLabel } from './plugins/rehypeTableLabel'

// 创建统一的 markdown 处理器
const processor = unified()
  .use(remarkParse)
  .use(remarkQQMusic)
  .use(remarkMath)
  .use(remarkGfm)
  .use(remarkRehype, {
    allowDangerousHtml: true
  })
  .use(rehypePrism, {
    showLineNumbers: true,
    ignoreMissing: true,
  })
  .use(rehypeTableLabel)
  .use(rehypeImgSize, {
    dir: path.join(process.cwd(), 'public')
  })
  .use(rehypeKatex, {
    strict: false
  })
  .use(rehypeStringify, {
    allowDangerousHtml: true
  })

// 异步渲染 markdown
async function renderMarkdown(content: string): Promise<string> {
  const result = await processor.process(content)
  return result.toString()
}

// 接口定义保持不变
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
  // 如果是字符串（单行格式），按逗号或空格分割
  if (typeof tags === 'string') {
    return tags
      .split(/[,\s]+/)
      .map(tag => tag.trim())
      .filter(Boolean)
  }
  
  // 如果是数组（YAML 列表格式），确保每个元素都是字符串
  if (Array.isArray(tags)) {
    return tags
      .map(tag => String(tag).trim())
      .filter(Boolean)
  }
  
  // 其他情况返回空数组
  return []
}

// 修改为异步函数
export async function getAllPosts(): Promise<BlogPost[]> {
  return getAllPostsFromDirectory('content/posts')
}

// 通用的从指定目录获取所有文章的函数
export async function getAllPostsFromDirectory(directory: string): Promise<BlogPost[]> {
  const postsDirectory = path.join(process.cwd(), directory)

  // 检查目录是否存在
  if (!fs.existsSync(postsDirectory)) {
    return []
  }

  const fileNames = fs.readdirSync(postsDirectory)

  const posts = await Promise.all(fileNames
    .filter(fileName => fileName.endsWith('.md'))
    .map(async fileName => {
      const { date, slug } = parseFileName(fileName)
      const fullPath = path.join(postsDirectory, fileName)
      const fileContents = fs.readFileSync(fullPath, 'utf8')
      const { data, content } = matter(fileContents)

      return {
        slug,
        title: data.title || slug,
        date,
        tags: parseTags(data.tags),
        excerpt: generateExcerpt(content),
        content: await renderMarkdown(content)
      }
    }))

  return posts.sort((a, b) => b.date.localeCompare(a.date))
}

// 获取 Schema 进度文章
export async function getAllSchemaProgress(): Promise<BlogPost[]> {
  return getAllPostsFromDirectory('content/schema-progress')
}

// 分页相关接口定义保持不变
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

// 修改为异步函数
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