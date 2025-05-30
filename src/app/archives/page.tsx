import { ArchivesList } from '@/components/archives/ArchivesList'
import { getAllPosts } from '@/lib/markdown'

export default async function ArchivesPage() {
  const posts = await getAllPosts()
  
  return (
    <div className="max-w-4xl mx-auto px-4 py-12">
      <h1 className="text-3xl font-bold mb-8 text-gray-900 dark:text-gray-100">
        归档
      </h1>
      <ArchivesList initialPosts={posts} />
    </div>
  )
} 