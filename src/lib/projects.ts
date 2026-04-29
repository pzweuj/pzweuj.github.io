import fs from 'fs'
import path from 'path'
import matter from 'gray-matter'
import { unified } from 'unified'
import remarkParse from 'remark-parse'
import remarkMath from 'remark-math'
import remarkRehype from 'remark-rehype'
import rehypeKatex from 'rehype-katex'
import rehypeStringify from 'rehype-stringify'
import rehypePrism from 'rehype-prism-plus'
import rehypeImgSize from 'rehype-img-size'
import remarkGfm from 'remark-gfm'
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
    aliases: {
      typescript: ['ts'],
      javascript: ['js'],
      python: ['py'],
      r: ['R'],
      perl: ['pl']
    }
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

let projectDocsCache: ProjectChapter[] | null = null

export async function getProjectDocs(): Promise<ProjectChapter[]> {
  if (projectDocsCache) return projectDocsCache

  const projectPath = path.join(process.cwd(), 'content/project')
  const chapters: ProjectChapter[] = []
  const htmlCache = loadHtmlCache()

  const dirs = fs.readdirSync(projectPath)

  for (const dir of dirs) {
    if (!dir.match(/^C\d+_/)) continue

    const chapterPath = path.join(projectPath, dir)
    const stat = fs.statSync(chapterPath)

    if (!stat.isDirectory()) continue

    const chapterOrder = parseInt(dir.split('_')[0].substring(1))
    const chapterTitle = dir.split('_')[1]

    const docs: ProjectDoc[] = []
    const files = fs.readdirSync(chapterPath)

    for (const file of files) {
      if (!file.endsWith('.md')) continue

      const docPath = path.join(chapterPath, file)
      const cacheKey = path.relative(process.cwd(), docPath).replace(/\\/g, '/')
      const content = fs.readFileSync(docPath, 'utf-8')
      const { data, content: markdown } = matter(content)

      const h1Title = extractH1Title(markdown)
      const contentWithoutTitle = markdown.replace(/^#\s+(.+)$/m, '')

      docs.push({
        slug: `${dir}/${file.replace('.md', '')}`,
        title: h1Title || data.title || file.replace('.md', '').split('_')[1] || file.replace('.md', ''),
        content: await renderMarkdownWithCache(contentWithoutTitle, cacheKey, 'project', htmlCache, renderMarkdown),
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
  saveHtmlCache(htmlCache)
  projectDocsCache = chapters
  return chapters
}

// 获取首页内容
let projectIndexCache: { title: string; content: string } | null = null

export async function getProjectIndex() {
  if (projectIndexCache) return projectIndexCache

  const indexPath = path.join(process.cwd(), 'content/project/index.md')
  const htmlCache = loadHtmlCache()
  const cacheKey = 'content/project/index.md'
  const content = fs.readFileSync(indexPath, 'utf-8')
  const { data, content: markdown } = matter(content)

  const result = {
    title: data.title || '实践项目',
    content: await renderMarkdownWithCache(markdown, cacheKey, 'project', htmlCache, renderMarkdown)
  }
  saveHtmlCache(htmlCache)
  projectIndexCache = result
  return result
} 