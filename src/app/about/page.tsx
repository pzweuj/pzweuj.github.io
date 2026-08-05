import { getAboutContent } from '@/lib/about'
import MarkdownContent from '@/components/blog/MarkdownContent'

export default async function AboutPage() {
  const about = await getAboutContent()

  return (
    <div className="max-w-4xl mx-auto px-4 py-12">
      <MarkdownContent html={about.content} />
    </div>
  )
} 