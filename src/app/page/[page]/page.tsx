import { getPaginatedPosts, getTotalPages } from '@/lib/markdown'
import { PostCard, Pagination } from '@/components/blog/PostList'
import { notFound } from 'next/navigation'

// 添加 generateStaticParams 函数
export async function generateStaticParams() {
  const totalPages = await getTotalPages() // 假设这个函数返回总页数
  
  return Array.from({ length: totalPages }, (_, i) => ({
    page: String(i + 1) // 将页码转换为字符串
  }))
}

export default async function BlogPage({ params }: { params: { page: string } }) {
  const { page: pageParam } = params
  const page = parseInt(pageParam)
  
  if (isNaN(page) || page < 1) {
    notFound()
  }

  const { posts, pagination } = await getPaginatedPosts(page)
  
  if (posts.length === 0) {
    notFound()
  }

  return (
    <div className="max-w-4xl mx-auto px-4 py-12">
      <div className="space-y-8">
        {/* 页面标题 */}
        {/* <div className="text-center space-y-2">
          <h1 className="text-3xl font-bold text-gray-900 dark:text-gray-100">博客文章</h1>
          <p className="text-gray-600 dark:text-gray-400">第 {page} 页</p>
        </div> */}

        {/* 文章列表 */}
        <div className="divide-y divide-gray-200 dark:divide-gray-800">
          {posts.map(post => (
            <PostCard key={post.slug} post={post} />
          ))}
        </div>

        {/* 分页 */}
        <Pagination {...pagination} />
      </div>
    </div>
  )
} 