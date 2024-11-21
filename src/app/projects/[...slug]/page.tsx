import { getProjectDocs } from '@/lib/projects'
import { ProjectSidebar } from '@/components/projects/ProjectSidebar'
import { notFound } from 'next/navigation'
import { Metadata } from 'next'

interface Props {
  params: {
    slug: string[]
  }
}

export default async function ProjectDocPage({ params }: Props) {
  const chapters = await getProjectDocs()
  const slug = params.slug.join('/')
  
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
        <div className="prose dark:prose-invert max-w-none">
          <h1>{doc.title}</h1>
          <div dangerouslySetInnerHTML={{ __html: doc.content }} />
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
  const chapters = await getProjectDocs()
  const slug = params.slug.join('/')
  
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