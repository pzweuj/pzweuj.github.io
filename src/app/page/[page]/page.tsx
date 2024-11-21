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
export default function BlogPage({ params }: PageProps) {
  const page = parseInt(params.page) // 将参数转为数字
  
  if (isNaN(page) || page < 1) {
    notFound() // 无效页码返回 404
  }

  const { posts, pagination } = getPaginatedPosts(page) // 获取分页数据
  
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
export function generateStaticParams() {
  const { pagination } = getPaginatedPosts(1) // 获取总页数
  const paths = []

  for (let i = 2; i <= pagination.totalPages; i++) {
    paths.push({ params: { page: i.toString() } }) // 包装为 params 格式
  }

  return paths
}
