import { unified } from 'unified'
import remarkParse from 'remark-parse'
import remarkMath from 'remark-math'
import remarkGfm from 'remark-gfm'
import remarkRehype from 'remark-rehype'
import rehypeKatex from 'rehype-katex'
import rehypeStringify from 'rehype-stringify'
import rehypeRaw from 'rehype-raw'
import rehypeSanitize, { defaultSchema } from 'rehype-sanitize'
import rehypePrism from 'rehype-prism-plus'
import { rehypeQQMusic } from './plugins/rehypeQQMusic'
import { rehypeTableLabel } from './plugins/rehypeTableLabel'
import { rehypeCodeCopy } from './plugins/rehypeCodeCopy'
import { rehypeLazyImage } from './plugins/rehypeLazyImage'
import { rehypeMermaid } from './plugins/rehypeMermaid'

// 插件链版本号，插件配置变更时递增以使缓存失效
export const PROCESSOR_VERSION = 4

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
    .use(remarkMath)
    .use(remarkGfm)
    .use(remarkRehype, { allowDangerousHtml: true })
    .use(rehypeRaw)
    .use(rehypeSanitize, {
      ...defaultSchema,
      tagNames: [...(defaultSchema.tagNames ?? []), 'font'],
      protocols: {
        ...defaultSchema.protocols,
        href: [...(defaultSchema.protocols?.href ?? []), 'ftp'],
      },
      // remark-rehype already namespaces footnotes. Keep their fragment links intact.
      clobber: ['name'],
      attributes: {
        ...defaultSchema.attributes,
        '*': [
          ...(defaultSchema.attributes?.['*'] ?? []).filter(attribute => attribute !== 'id'),
          ['id', /^user-content-/, 'footnote-label'],
        ],
        font: [['color', 'red']],
        code: [['className', /^language-./, 'math-inline', 'math-display']],
      },
    })
    .use(rehypeQQMusic)
    .use(rehypeMermaid)
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
    .use(rehypeStringify)
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
    return '<p role="alert">这篇文章暂时无法渲染。</p>'
  }
}
