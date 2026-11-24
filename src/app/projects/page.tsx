import { getProjectIndex } from '@/lib/projects'
import MarkdownContent from '@/components/blog/MarkdownContent'

export default async function ProjectsPage() {
  const index = await getProjectIndex()
  return <MarkdownContent html={index.content} />
}
