import { redirect } from 'next/navigation'

import { getPaginatedPosts } from '@/lib/markdown'
import { PostCard, Pagination } from '@/components/blog/PostList'

export default async function Home() {
  const { posts, pagination } = await getPaginatedPosts(1)

  return (
    <div className="max-w-4xl mx-auto px-4 py-12">
      <div className="space-y-8">
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
