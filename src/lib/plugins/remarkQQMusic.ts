import { visit } from 'unist-util-visit'
import type { Plugin } from 'unified'
import type { Root } from 'mdast'

interface QQMusicOptions {
  className?: string
}

const QQ_MUSIC_REGEX = /^\[qqmusic:(\d+)\]$/

export const remarkQQMusic: Plugin<[QQMusicOptions?], Root> = (options = {}) => {
  const className = options.className || 'qq-music-player'
  
  return (tree) => {
    visit(tree, 'paragraph', (node: any) => {
      if (node.children.length === 1 && node.children[0].type === 'text') {
        const match = node.children[0].value.match(QQ_MUSIC_REGEX)
        if (match) {
          const songId = match[1]
          node.type = 'html'
          node.children = undefined
          node.value = `<div class="${className}"><iframe frameborder="no" border="0" marginwidth="0" marginheight="0" width="100%" height="86" src="https://i.y.qq.com/n2/m/outchain/player/index.html?songid=${songId}"></iframe></div>`
        }
      }
    })
  }
} 