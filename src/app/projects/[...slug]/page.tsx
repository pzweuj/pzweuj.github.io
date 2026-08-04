import { getProjectDocs } from '@/lib/projects'
import { ProjectSidebar } from '@/components/projects/ProjectSidebar'
import { notFound } from 'next/navigation'
import { Metadata } from 'next'
import MarkdownContent from '@/components/blog/MarkdownContent'

interface Props {
  params: Promise<{
    slug: string[]
  }>
}

export default async function ProjectDocPage({ params }: Props) {
  const { slug: slugArr } = await params
  const chapters = await getProjectDocs()
  const slug = slugArr.join('/')
  
  // 查找当前文档
  const doc = chapters
    .flatMap(chapter => chapter.docs)
    .find(doc => doc.slug === slug)
  
  if (!doc) {
    notFound()
  }
  
  return (
    <div className="flex">
      <ProjectSidebar chapters={chapters} />
      
      <main className="flex-1 max-w-4xl mx-auto px-4 py-12">
        <div className="prose max-w-none">
          <h1>{doc.title}</h1>
          <MarkdownContent html={doc.content} className="prose max-w-none" />
        </div>
      </main>
    </div>
  )
}

export async function generateStaticParams() {
  const chapters = await getProjectDocs()
  return chapters.flatMap(chapter => 
    chapter.docs.map(doc => ({
      slug: doc.slug.split('/')
    }))
  )
}

export async function generateMetadata({ params }: Props): Promise<Metadata> {
  const { slug: slugArr } = await params
  const chapters = await getProjectDocs()
  const slug = slugArr.join('/')

  const doc = chapters
    .flatMap(chapter => chapter.docs)
    .find(doc => doc.slug === slug)
  
  if (!doc) {
    return {
      title: '文档未找到',
      description: '请求的项目文档不存在'
    }
  }

  return {
    title: doc.title,
    description: `${doc.title} - 项目文档`
  }
} 