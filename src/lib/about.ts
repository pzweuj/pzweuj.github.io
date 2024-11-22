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

// 创建统一的 markdown 处理器
const processor = unified()
  .use(remarkParse)
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

export async function getAboutContent() {
  const filePath = path.join(process.cwd(), 'content/about.md')
  const content = fs.readFileSync(filePath, 'utf-8')
  const { data, content: markdown } = matter(content)
  
  return {
    title: data.title || '关于',
    content: await renderMarkdown(markdown)
  }
} 