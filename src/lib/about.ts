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
import { loadHtmlCache, saveHtmlCache, renderMarkdownWithCache } from './cache'

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
    content: await renderMarkdownWithCache(markdown, cacheKey, 'about', htmlCache, renderMarkdown)
  }
  saveHtmlCache(htmlCache)
  aboutCache = result
  return result
} 