import type { ReactNode } from 'react'
import { ProjectSidebar } from '@/components/projects/ProjectSidebar'
import { getProjectDocs } from '@/lib/projects'

export default async function ProjectsLayout({ children }: { children: ReactNode }) {
  const chapters = (await getProjectDocs()).map(({ id, title, docs }) => ({
    id,
    title,
    docs: docs.map(({ slug, title }) => ({ slug, title })),
  }))

  return (
    <div className="mx-auto max-w-7xl md:flex">
      <ProjectSidebar chapters={chapters} />
      <div className="mx-auto min-w-0 max-w-4xl flex-1 px-4 py-8 sm:px-6 md:px-8 md:py-12">
        {children}
      </div>
    </div>
  )
}
