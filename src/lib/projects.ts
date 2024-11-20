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
import 'prismjs/components/prism-python'  // Python 支持
import 'prismjs/components/prism-r'       // R 语言支持
import 'prismjs/components/prism-perl'    // Perl 支持

const md = new MarkdownIt({
  html: true,
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

export interface ProjectDoc {
  slug: string
  title: string
  content: string
  chapter: string
  order: number
}

export interface ProjectChapter {
  id: string
  title: string
  order: number
  docs: ProjectDoc[]
}

// 从 Markdown 内容中提取第一个 h1 标题
function extractH1Title(content: string): string | null {
  const match = content.match(/^#\s+(.+)$/m)
  return match ? match[1].trim() : null
}

export async function getProjectDocs() {
  const projectPath = path.join(process.cwd(), 'content/project')
  const chapters: ProjectChapter[] = []
  
  const dirs = await fs.readdir(projectPath)
  
  for (const dir of dirs) {
    if (!dir.match(/^C\d+_/)) continue
    
    const chapterPath = path.join(projectPath, dir)
    const stat = await fs.stat(chapterPath)
    
    if (!stat.isDirectory()) continue
    
    const chapterOrder = parseInt(dir.split('_')[0].substring(1))
    const chapterTitle = dir.split('_')[1]
    
    const docs: ProjectDoc[] = []
    const files = await fs.readdir(chapterPath)
    
    for (const file of files) {
      if (!file.endsWith('.md')) continue
      
      const docPath = path.join(chapterPath, file)
      const content = await fs.readFile(docPath, 'utf-8')
      const { data, content: markdown } = matter(content)
      
      // 优先使用 h1 标题，如果没有则使用文件名
      const h1Title = extractH1Title(markdown)
      // 移除 markdown 中的 h1 标题
      const contentWithoutTitle = markdown.replace(/^#\s+(.+)$/m, '')
      
      docs.push({
        slug: `${dir}/${file.replace('.md', '')}`,
        title: h1Title || data.title || file.replace('.md', '').split('_')[1] || file.replace('.md', ''),
        content: md.render(contentWithoutTitle), // 使用去除标题后的内容
        chapter: dir,
        order: parseInt(file.split('_')[0])
      })
    }
    
    docs.sort((a, b) => a.order - b.order)
    
    chapters.push({
      id: dir,
      title: chapterTitle,
      order: chapterOrder,
      docs
    })
  }
  
  chapters.sort((a, b) => a.order - b.order)
  
  return chapters
}

// 获取首页内容
export async function getProjectIndex() {
  const indexPath = path.join(process.cwd(), 'content/project/index.md')
  const content = await fs.readFile(indexPath, 'utf-8')
  const { data, content: markdown } = matter(content)
  
  return {
    title: data.title || '实践项目',
    content: md.render(markdown)
  }
} 