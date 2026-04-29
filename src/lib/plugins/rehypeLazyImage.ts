import { visit } from 'unist-util-visit'
import type { Element, Parent } from 'hast'

/**
 * rehype 插件：为所有 img 标签添加 loading="lazy" 和 decoding="async"
 */
export function rehypeLazyImage() {
  return (tree: Parent) => {
    visit(tree, 'element', (node: Element) => {
      if (node.tagName !== 'img') return
      node.properties = node.properties || {}
      if (!node.properties.loading) {
        node.properties.loading = 'lazy'
      }
      if (!node.properties.decoding) {
        node.properties.decoding = 'async'
      }
    })
  }
}
