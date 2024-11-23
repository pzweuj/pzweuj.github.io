import { visit } from 'unist-util-visit'
import type { Plugin } from 'unified'
import type { Root, Paragraph, Text } from 'mdast'

interface QQMusicOptions {
  className?: string
}

const QQ_MUSIC_REGEX = /^\[qqmusic:(\d+)\]$/

export const remarkQQMusic: Plugin<[QQMusicOptions?], Root> = (options = {}) => {
  const className = options.className || 'qq-music-player'
  
  return (tree) => {
    visit(tree, 'paragraph', (node: Paragraph) => {
      if (
        node.children.length === 1 && 
        node.children[0].type === 'text'
      ) {
        const textNode = node.children[0] as Text
        const match = textNode.value.match(QQ_MUSIC_REGEX)
        if (match) {
          const songId = match[1]
          const htmlNode = node as unknown as { 
            type: 'html'
            value: string
            children?: never
          }
          htmlNode.type = 'html'
          htmlNode.children = undefined
          htmlNode.value = `<div class="${className}"><iframe frameborder="no" border="0" marginwidth="0" marginheight="0" width="100%" height="86" src="https://i.y.qq.com/n2/m/outchain/player/index.html?songid=${songId}"></iframe></div>`
        }
      }
    })
  }
} 