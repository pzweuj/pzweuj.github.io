declare module 'rehype-prism-plus' {
  import type { Plugin } from 'unified'

  interface RehypePrismOptions {
    showLineNumbers?: boolean
    ignoreMissing?: boolean
    aliases?: {
      [key: string]: string[]
    }
  }

  const rehypePrism: Plugin<[RehypePrismOptions?]>
  export default rehypePrism
} 