import { BlogPost, getAllPosts } from './markdown'

// 按年份组织文章
export function getPostsByYear() {
  const posts = getAllPosts()
  const postsByYear: { [key: string]: BlogPost[] } = {}

  posts.forEach(post => {
    const year = post.date.substring(0, 4)
    if (!postsByYear[year]) {
      postsByYear[year] = []
    }
    postsByYear[year].push(post)
  })

  // 返回按年份排序的对象
  return Object.entries(postsByYear)
    .sort(([yearA], [yearB]) => yearB.localeCompare(yearA))
    .reduce((acc, [year, posts]) => {
      acc[year] = posts
      return acc
    }, {} as { [key: string]: BlogPost[] })
}

// 获取所有标签
export function getAllTags() {
  const posts = getAllPosts()
  const tags = new Set<string>()
  
  posts.forEach(post => {
    post.tags.forEach(tag => tags.add(tag))
  })

  return Array.from(tags).sort()
}

// 按标签筛选文章
export function getPostsByTag(tag: string) {
  const posts = getAllPosts()
  return posts.filter(post => post.tags.includes(tag))
} 