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
import remarkDirective from 'remark-directive'
import { visit } from 'unist-util-visit'

// 添加 QQ 音乐解析插件
function remarkQQMusic() {
  return (tree: any) => {
    visit(tree, (node) => {
      if (
        node.type === 'textDirective' ||
        node.type === 'leafDirective' ||
        node.type === 'containerDirective'
      ) {
        if (node.name !== 'qqmusic') return

        const data = node.data || (node.data = {})
        const songId = node.attributes?.id

        if (!songId) return

        data.hName = 'div'
        data.hProperties = {
          className: 'qq-music-player',
          'data-song-id': songId,
        }
        
        const iframeHtml = `<iframe 
          frameborder="no" 
          border="0" 
          marginwidth="0" 
          marginheight="0" 
          width="100%" 
          height="86" 
          src="https://i.y.qq.com/n2/m/outchain/player/index.html?songid=${songId}&songtype=0"
        ></iframe>`

        data.hChildren = [{
          type: 'raw',
          value: iframeHtml
        }]
      }
    })
  }
}

// 创建统一的 markdown 处理器
const processor = unified()
  .use(remarkParse)
  .use(remarkDirective)
  .use(remarkQQMusic)
  .use(remarkMath)
  .use(remarkRehype, { allowDangerousHtml: true })
  .use(rehypePrism, {
    showLineNumbers: true,
    ignoreMissing: true,
  })
  .use(rehypeKatex, {
    strict: false
  })
  .use(rehypeStringify, { allowDangerousHtml: true })

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