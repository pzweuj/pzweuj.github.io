'use client'

import { useState, useSyncExternalStore } from 'react'
import Link from 'next/link'
import { usePathname } from 'next/navigation'
import { BookOpen, ChevronDown, ListTree, PanelLeftClose, PanelLeftOpen } from 'lucide-react'
import type { ProjectChapter, ProjectDoc } from '@/lib/projects'

type Chapter = Pick<ProjectChapter, 'id' | 'title'> & {
  docs: Pick<ProjectDoc, 'slug' | 'title'>[]
}

function sidebarPreference(key: string) {
  let unsaved: string | null = null
  const listeners = new Set<() => void>()

  return {
    read() {
      if (unsaved !== null) return unsaved
      try {
        return localStorage.getItem(key)
      } catch {
        return null
      }
    },
    write(value: string) {
      try {
        localStorage.setItem(key, value)
        unsaved = null
      } catch {
        unsaved = value
      }
      listeners.forEach(notify => notify())
    },
    subscribe(notify: () => void) {
      const onStorage = (event: StorageEvent) => {
        if (event.key === key || event.key === null) {
          unsaved = null
          notify()
        }
      }
      listeners.add(notify)
      window.addEventListener('storage', onStorage)
      return () => {
        listeners.delete(notify)
        window.removeEventListener('storage', onStorage)
      }
    },
  }
}

const expanded = sidebarPreference('expandedChapters')
const collapsed = sidebarPreference('projectSidebarCollapsed')
const serverSnapshot = () => null

function readChapterIds(value: string | null): string[] {
  try {
    const parsed: unknown = JSON.parse(value ?? '[]')
    return Array.isArray(parsed) ? parsed.filter(id => typeof id === 'string') : []
  } catch {
    return []
  }
}

const activeLink = 'bg-blue-50 font-medium text-blue-700 dark:bg-blue-950 dark:text-blue-300'
const idleLink = 'text-gray-600 hover:bg-gray-100 hover:text-gray-900 dark:text-gray-300 dark:hover:bg-gray-800 dark:hover:text-white'

export function ProjectSidebar({ chapters }: { chapters: Chapter[] }) {
  const pathname = usePathname().replace(/\/$/, '').replace(/\.html$/, '')
  return <Directory key={pathname} pathname={pathname} chapters={chapters} />
}

function Directory({ pathname, chapters }: { pathname: string; chapters: Chapter[] }) {
  const isCollapsed = useSyncExternalStore(collapsed.subscribe, collapsed.read, serverSnapshot) === 'true'
  const saved = useSyncExternalStore(expanded.subscribe, expanded.read, serverSnapshot)
  const expandedIds = readChapterIds(saved)
  const [mobileOpen, setMobileOpen] = useState(false)
  const [activeChapterClosed, setActiveChapterClosed] = useState(false)
  const activeChapter = pathname.split('/')[2]
  const toggleLabel = isCollapsed ? '展开目录' : '收起目录'
  const DesktopIcon = isCollapsed ? PanelLeftOpen : PanelLeftClose

  function toggleChapter(id: string, open: boolean) {
    if (id === activeChapter) setActiveChapterClosed(open)
    const next = expandedIds.filter(chapterId => chapterId !== id)
    if (!open) next.push(id)
    expanded.write(JSON.stringify(next))
  }

  return (
    <aside className={`sticky top-14 z-30 flex flex-col self-start border-b border-gray-200 bg-white dark:border-gray-800 dark:bg-gray-900 md:top-16 md:h-[calc(100dvh-4rem)] md:shrink-0 md:border-r md:border-b-0 ${isCollapsed ? 'md:w-14' : 'md:w-64'}`}>
      <button
        type="button"
        aria-expanded={mobileOpen}
        aria-controls="project-directory"
        onClick={() => setMobileOpen(!mobileOpen)}
        className="flex min-h-11 items-center gap-2 px-4 text-sm font-semibold text-gray-800 dark:text-gray-100 md:hidden"
      >
        <ListTree size={18} aria-hidden />
        实践目录
        <ChevronDown size={16} aria-hidden className={`ml-auto ${mobileOpen ? 'rotate-180' : ''}`} />
      </button>
      <div className={`hidden min-h-14 items-center md:flex ${isCollapsed ? 'justify-center' : 'justify-between px-3'}`}>
        {!isCollapsed && (
          <span className="flex items-center gap-2 pl-2 text-sm font-semibold text-gray-900 dark:text-gray-100">
            <ListTree size={18} aria-hidden /> 实践目录
          </span>
        )}
        <button
          type="button"
          aria-label={toggleLabel}
          title={toggleLabel}
          aria-expanded={!isCollapsed}
          aria-controls="project-directory"
          onClick={() => collapsed.write(String(!isCollapsed))}
          className={`grid size-9 place-items-center rounded-md ${idleLink}`}
        >
          <DesktopIcon size={18} aria-hidden />
        </button>
      </div>
      <nav
        id="project-directory"
        aria-label="实践目录"
        className={`max-h-[calc(100dvh-7rem)] overflow-y-auto overscroll-contain px-3 pb-4 text-sm md:max-h-none ${mobileOpen ? 'block' : 'hidden'} ${isCollapsed ? 'md:hidden' : 'md:block'}`}
      >
        <Link
          href="/projects"
          aria-current={pathname === '/projects' ? 'page' : undefined}
          onClick={() => setMobileOpen(false)}
          className={`flex min-h-10 items-center gap-2 rounded-md px-2 ${pathname === '/projects' ? activeLink : idleLink}`}
        >
          <BookOpen size={16} aria-hidden /> 实践首页
        </Link>
        <ul className="mt-2 space-y-1 border-t border-gray-200 pt-2 dark:border-gray-800">
          {chapters.map(chapter => {
            const current = chapter.id === activeChapter
            const open = current ? !activeChapterClosed : expandedIds.includes(chapter.id)
            const listId = `project-docs-${chapter.id}`
            return (
              <li key={chapter.id}>
                <button
                  type="button"
                  aria-expanded={open}
                  aria-controls={listId}
                  onClick={() => toggleChapter(chapter.id, open)}
                  className="flex min-h-11 w-full items-center gap-2 rounded-md px-2 text-left font-semibold text-gray-800 hover:bg-gray-100 dark:text-gray-200 dark:hover:bg-gray-800 md:min-h-9"
                >
                  <ChevronDown size={16} aria-hidden className={`shrink-0 text-gray-500 motion-safe:transition-transform ${open ? '' : '-rotate-90'}`} />
                  <span className="min-w-0 break-words">{chapter.title}</span>
                </button>
                <ul id={listId} hidden={!open} className="my-1 ml-4 space-y-0.5 border-l border-gray-200 pl-2 dark:border-gray-700">
                  {chapter.docs.map(doc => {
                    const href = `/projects/${doc.slug}`
                    const selected = pathname === href
                    return (
                      <li key={doc.slug}>
                        <Link
                          href={href}
                          aria-current={selected ? 'page' : undefined}
                          onClick={() => setMobileOpen(false)}
                          className={`block break-words rounded-md px-2.5 py-2 leading-5 md:py-1.5 ${selected ? activeLink : idleLink}`}
                        >
                          {doc.title}
                        </Link>
                      </li>
                    )
                  })}
                </ul>
              </li>
            )
          })}
        </ul>
      </nav>
    </aside>
  )
}
