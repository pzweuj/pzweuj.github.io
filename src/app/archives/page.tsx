import { ArchivesList } from '@/components/archives/ArchivesList'
import { getAllPosts } from '@/lib/markdown'

export default async function ArchivesPage() {
  const posts = await getAllPosts()
  
  return (
    <div className="max-w-4xl mx-auto px-4 py-12">
      {/* <div className="text-center mb-8">
        <h1 className="text-3xl font-bold text-gray-900 dark:text-gray-100">
          文章归档
        </h1>
        <p className="mt-2 text-gray-600 dark:text-gray-400">
          共 {posts.length} 篇文章
        </p>
      </div> */}
      
      <ArchivesList initialPosts={posts} />
    </div>
  )
} 