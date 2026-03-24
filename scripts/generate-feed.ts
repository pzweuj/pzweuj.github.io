import fs from 'fs'
import path from 'path'
import matter from 'gray-matter'
import { Feed } from 'feed'
import selfConfig from '../src/config/self.config'

interface Post {
  slug: string
  title: string
  date: string
  tags: string[]
  excerpt: string
  content: string
}

// 从文件名中提取日期和slug
function parseFileName(fileName: string) {
  const match = fileName.match(/^(\d{4}-\d{2}-\d{2})-(.+)\.md$/)
  if (!match) return null
  return {
    date: match[1],
    slug: match[2]
  }
}

// 生成摘要
function generateExcerpt(content: string, maxLength: number = 200) {
  const plainText = content
    .replace(/<[^>]+>/g, '')
    .replace(/[#*`]/g, '')
    .trim()

  const firstParagraph = plainText
    .split('\n')
    .find(line => line.trim().length > 0)
    ?.trim() || ''

  return firstParagraph.length > maxLength
    ? `${firstParagraph.slice(0, maxLength)}...`
    : firstParagraph
}

function parseTags(tags: unknown): string[] {
  if (typeof tags === 'string') {
    return tags
      .split(/[,\s]+/)
      .map(tag => tag.trim())
      .filter(Boolean)
  }

  if (Array.isArray(tags)) {
    return tags
      .map(tag => String(tag).trim())
      .filter(Boolean)
  }

  return []
}

// 获取所有文章（不渲染 markdown，只获取元数据）
function getAllPosts(): Post[] {
  const postsDirectory = path.join(process.cwd(), 'content/posts')

  if (!fs.existsSync(postsDirectory)) {
    return []
  }

  const fileNames = fs.readdirSync(postsDirectory)

  const posts = fileNames
    .filter(fileName => fileName.endsWith('.md'))
    .map(fileName => {
      const parsed = parseFileName(fileName)
      if (!parsed) return null

      const { date, slug } = parsed
      const fullPath = path.join(postsDirectory, fileName)
      const fileContents = fs.readFileSync(fullPath, 'utf8')
      const { data, content } = matter(fileContents)

      return {
        slug,
        title: data.title || slug,
        date,
        tags: parseTags(data.tags),
        excerpt: generateExcerpt(content),
        content: content // 保留原始内容用于 feed
      }
    })
    .filter((post): post is Post => post !== null)

  return posts.sort((a, b) => b.date.localeCompare(a.date))
}

async function generateFeed() {
  const posts = getAllPosts()
  const siteUrl = selfConfig.siteUrl

  const feed = new Feed({
    title: selfConfig.title,
    description: selfConfig.description,
    id: siteUrl,
    link: siteUrl,
    language: selfConfig.language,
    image: `${siteUrl}/icon.ico`,
    favicon: `${siteUrl}/icon.ico`,
    copyright: `Copyright © ${new Date().getFullYear()} ${selfConfig.author}`,
    updated: new Date(),
    feedLinks: {
      rss2: `${siteUrl}/rss.xml`,
      atom: `${siteUrl}/atom.xml`,
      json: `${siteUrl}/feed.json`,
    },
    author: {
      name: selfConfig.author,
      email: selfConfig.social.email,
      link: selfConfig.social.github,
    },
  })

  posts.forEach(post => {
    const url = `${siteUrl}/posts/${post.slug}`

    feed.addItem({
      title: post.title,
      id: url,
      link: url,
      description: post.excerpt,
      content: post.excerpt,
      date: new Date(post.date),
      category: post.tags.map(tag => ({ name: tag })),
    })
  })

  // 写入文件到 public 目录
  const publicDir = path.join(process.cwd(), 'public')

  fs.writeFileSync(path.join(publicDir, 'rss.xml'), feed.rss2())
  fs.writeFileSync(path.join(publicDir, 'atom.xml'), feed.atom1())
  fs.writeFileSync(path.join(publicDir, 'feed.json'), feed.json1())

  console.log('Feed files generated successfully!')
  console.log(`- rss.xml (${posts.length} posts)`)
  console.log(`- atom.xml (${posts.length} posts)`)
  console.log(`- feed.json (${posts.length} posts)`)
}

generateFeed().catch(console.error)