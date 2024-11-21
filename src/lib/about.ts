import fs from 'fs'
import path from 'path'
import matter from 'gray-matter'
import MarkdownIt from 'markdown-it'

const md = new MarkdownIt({
  html: true,
  breaks: true,
  linkify: true,
  typographer: true
})

export function getAboutContent() {
  const filePath = path.join(process.cwd(), 'content/about.md')
  const content = fs.readFileSync(filePath, 'utf-8')
  const { data, content: markdown } = matter(content)
  
  return {
    title: data.title || '关于',
    content: md.render(markdown)
  }
} 