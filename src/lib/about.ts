'use server'

import fs from 'fs/promises'
import path from 'path'
import matter from 'gray-matter'
import MarkdownIt from 'markdown-it'

const md = new MarkdownIt({
  html: true,
  breaks: true,
  linkify: true,
  typographer: true
})

export async function getAboutContent() {
  const filePath = path.join(process.cwd(), 'content/about.md')
  const content = await fs.readFile(filePath, 'utf-8')
  const { data, content: markdown } = matter(content)
  
  return {
    title: data.title || '关于',
    content: md.render(markdown)
  }
} 