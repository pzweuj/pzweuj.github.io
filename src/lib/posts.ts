import { getAllPosts, BlogPost } from '@/lib/markdown'

export interface PostsByYear {
  [key: string]: BlogPost[]
}

export async function getPostsByYear(): Promise<PostsByYear> {
  const posts = await getAllPosts()
  const postsByYear: PostsByYear = {}

  posts.forEach(post => {
    const year = post.date.substring(0, 4)
    if (!postsByYear[year]) {
      postsByYear[year] = []
    }
    postsByYear[year].push(post)
  })

  // 对每年的文章按日期排序
  Object.keys(postsByYear).forEach(year => {
    postsByYear[year].sort((a, b) => b.date.localeCompare(a.date))
  })

  return postsByYear
}

// 如果有其他函数也需要处理异步数据，同样添加 async/await
export async function getPostTags(): Promise<string[]> {
  const posts = await getAllPosts()
  const tagSet = new Set<string>()
  
  posts.forEach(post => {
    post.tags.forEach(tag => tagSet.add(tag))
  })

  return Array.from(tagSet).sort()
}

// 其他相关函数也需要类似处理... 