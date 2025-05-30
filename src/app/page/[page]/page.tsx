import { getPaginatedPosts } from '@/lib/markdown'
import { PostCard, Pagination } from '@/components/blog/PostList'
import { notFound } from 'next/navigation'

// 定义类型
type PageProps = {
  params: {
    page: string
  }
  searchParams?: { [key: string]: string | string[] | undefined }
}

// 页面组件
export default async function BlogPage({ params }: PageProps) {
  const page = parseInt(params.page) // 将参数转为数字
  
  if (isNaN(page) || page < 1) {
    notFound() // 无效页码返回 404
  }

  const { posts, pagination } = await getPaginatedPosts(page) // 添加 await
  
  if (posts.length === 0) {
    notFound() // 无文章内容返回 404
  }

  return (
    <div className="max-w-4xl mx-auto px-4 py-12">
      <div className="space-y-8">
        <div className="divide-y divide-gray-200 dark:divide-gray-800">
          {posts.map(post => (
            <PostCard key={post.slug} post={post} />
          ))}
        </div>
        <Pagination {...pagination} />
      </div>
    </div>
  )
}

// 生成静态路径
export async function generateStaticParams() {  // 添加 async
  const { pagination } = await getPaginatedPosts(1) // 添加 await
  
  // 生成所有可能的页码参数
  return Array.from({ length: pagination.totalPages }, (_, i) => ({
    page: (i + 1).toString()
  }))
}
