import { getAboutContent } from '@/lib/about'

export default async function AboutPage() {
  const about = await getAboutContent()
  
  return (
    <div className="max-w-4xl mx-auto px-4 py-12">
      <div className="prose dark:prose-invert max-w-none">
        <div dangerouslySetInnerHTML={{ __html: about.content }} />
      </div>
    </div>
  )
} 