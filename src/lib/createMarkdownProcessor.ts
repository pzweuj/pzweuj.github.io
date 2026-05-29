import { unified } from 'unified'
import remarkParse from 'remark-parse'
import remarkMath from 'remark-math'
import remarkGfm from 'remark-gfm'
import remarkRehype from 'remark-rehype'
import rehypeKatex from 'rehype-katex'
import rehypeStringify from 'rehype-stringify'
import rehypePrism from 'rehype-prism-plus'
import { remarkQQMusic } from './plugins/remarkQQMusic'
import { rehypeTableLabel } from './plugins/rehypeTableLabel'
import { rehypeCodeCopy } from './plugins/rehypeCodeCopy'
import { rehypeLazyImage } from './plugins/rehypeLazyImage'

// 插件链版本号，插件配置变更时递增以使缓存失效
export const PROCESSOR_VERSION = 1

// Prism.js 语言别名，统一用于所有管道
const PRISM_ALIASES: Record<string, string[]> = {
  typescript: ['ts'],
  javascript: ['js'],
  python: ['py'],
  r: ['R'],
  perl: ['pl'],
}

export interface ProcessorOptions {
  enableTableLabel?: boolean
  enableCodeCopy?: boolean
  enableLazyImage?: boolean
}

/**
 * 创建统一的 Markdown 处理器。
 * 工厂函数，每次调用创建新实例，避免模块顶层副作用。
 */
export function createMarkdownProcessor(options: ProcessorOptions = {}) {
  const {
    enableTableLabel = false,
    enableCodeCopy = false,
    enableLazyImage = false,
  } = options

  // 构建完整管道：remark 阶段 → rehype 阶段
  // 条件插件通过中间变量承接，避免 let 重赋值的类型推断问题
  const base = unified()
    .use(remarkParse)
    .use(remarkQQMusic)
    .use(remarkMath)
    .use(remarkGfm)
    .use(remarkRehype, { allowDangerousHtml: true })
    .use(rehypePrism, {
      showLineNumbers: true,
      ignoreMissing: true,
      aliases: PRISM_ALIASES,
    })

  // eslint-disable-next-line @typescript-eslint/no-explicit-any
  let proc: any = base

  if (enableTableLabel) proc = proc.use(rehypeTableLabel)
  if (enableCodeCopy) proc = proc.use(rehypeCodeCopy)
  if (enableLazyImage) proc = proc.use(rehypeLazyImage)

  return proc
    .use(rehypeKatex, { strict: false })
    .use(rehypeStringify, { allowDangerousHtml: true })
}

/**
 * 渲染 Markdown 为 HTML 字符串，带错误处理。
 * 单篇文章渲染失败不会导致全站构建崩溃。
 */
export async function renderMarkdown(
  processor: ReturnType<typeof createMarkdownProcessor>,
  content: string,
  filePath?: string,
): Promise<string> {
  try {
    const result = await processor.process(content)
    return result.toString()
  } catch (error) {
    const label = filePath ? `[${filePath}]` : '[unknown]'
    console.error(`Markdown 渲染失败 ${label}:`, error)
    return `<div style="color:red;border:1px solid red;padding:1em;margin:1em 0;">
      <strong>渲染失败</strong> ${label}<br/>
      <pre>${(error as Error).message}</pre>
    </div>`
  }
}
