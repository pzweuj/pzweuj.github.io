'use server'

import fs from 'fs/promises'
import path from 'path'
import matter from 'gray-matter'
import MarkdownIt from 'markdown-it'
import Prism from 'prismjs'
import katex from 'markdown-it-katex'

// 加载更多语言支持
import 'prismjs/components/prism-typescript'
import 'prismjs/components/prism-javascript'
import 'prismjs/components/prism-jsx'
import 'prismjs/components/prism-tsx'
import 'prismjs/components/prism-bash'
import 'prismjs/components/prism-json'
import 'prismjs/components/prism-markdown'
import 'prismjs/components/prism-css'
import 'prismjs/components/prism-python'
import 'prismjs/components/prism-r'
import 'prismjs/components/prism-perl'

const md = new MarkdownIt({
  html: true,
  breaks: true,
  linkify: true,
  typographer: true,
  highlight: (str, lang) => {
    if (lang && Prism.languages[lang]) {
      try {
        return `<pre class="language-${lang}"><code>${Prism.highlight(str, Prism.languages[lang], lang)}</code></pre>`
      } catch (err) {
        console.error(err)
      }
    }
    return `<pre class="language-${lang}"><code>${md.utils.escapeHtml(str)}</code></pre>`
  }
}).use(katex)  // 启用 KaTeX 支持

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

export async function getAllPosts(): Promise<BlogPost[]> {
  const postsDirectory = path.join(process.cwd(), 'content/posts')
  const fileNames = await fs.readdir(postsDirectory)
  
  const posts = await Promise.all(
    fileNames
      .filter(fileName => fileName.endsWith('.md'))
      .map(async fileName => {
        const { date, slug } = parseFileName(fileName)
        const fullPath = path.join(postsDirectory, fileName)
        const fileContents = await fs.readFile(fullPath, 'utf8')
        const { data, content } = matter(fileContents)
        
        return {
          slug,
          title: data.title || slug,
          date,
          tags: parseTags(data.tags),
          excerpt: generateExcerpt(content),
          content: md.render(content)
        }
      })
  )

  return posts.sort((a, b) => b.date.localeCompare(a.date))
}

export async function getPaginatedPosts(page: number = 1, limit: number = 5) {
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

// 假设每页显示 10 篇文章
const POSTS_PER_PAGE = 10

export async function getTotalPages() {
  const posts = await getAllPosts()
  return Math.ceil(posts.length / POSTS_PER_PAGE)
} 