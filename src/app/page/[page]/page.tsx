import { getPaginatedPosts } from '@/lib/markdown'
import { PostCard, Pagination } from '@/components/blog/PostList'
import { notFound } from 'next/navigation'

export default function BlogPage({ params }: { params: { page: string } }) {
  const page = parseInt(params.page)
  
  if (isNaN(page) || page < 1) {
    notFound()
  }

  const { posts, pagination } = getPaginatedPosts(page)
  
  if (posts.length === 0) {
    notFound()
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

export function generateStaticParams() {
  const { pagination } = getPaginatedPosts(1)
  const paths = []
  
  for (let i = 2; i <= pagination.totalPages; i++) {
    paths.push({ page: i.toString() })
  }
  
  return paths
} 