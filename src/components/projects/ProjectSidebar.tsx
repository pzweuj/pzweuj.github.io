'use client'

import { useState, useSyncExternalStore } from 'react'
import Link from 'next/link'
import { usePathname } from 'next/navigation'
import type { ProjectChapter } from '@/lib/projects'
import { Menu, X } from 'lucide-react'

const EXPANDED_KEY = 'expandedChapters'

// --- localStorage 外部 store（useSyncExternalStore 用）---
// 用模块级缓存保证 getSnapshot 在值未变化时返回同一引用，
// 否则 useSyncExternalStore 会无限重渲染。
let expandedCache: string[] | null = null
const expandedListeners = new Set<() => void>()

function readExpanded(): string[] {
  if (typeof window === 'undefined') return []
  try {
    const saved = localStorage.getItem(EXPANDED_KEY)
    if (!saved) return []
    const parsed = JSON.parse(saved)
    return Array.isArray(parsed) ? parsed : []
  } catch {
    return []
  }
}

function getExpandedSnapshot(): string[] {
  const next = readExpanded()
  if (expandedCache === null) {
    expandedCache = next
    return expandedCache
  }
  const changed =
    expandedCache.length !== next.length ||
    expandedCache.some((id, i) => id !== next[i])
  if (changed) expandedCache = next
  return expandedCache
}

function getExpandedServerSnapshot(): string[] {
  return []
}

function subscribeExpanded(callback: () => void) {
  expandedListeners.add(callback)
  const onStorage = (e: StorageEvent) => {
    if (e.key === EXPANDED_KEY) {
      expandedCache = null // 强制下次 getSnapshot 重读
      callback()
    }
  }
  window.addEventListener('storage', onStorage)
  return () => {
    expandedListeners.delete(callback)
    window.removeEventListener('storage', onStorage)
  }
}

function writeExpanded(next: string[]) {
  localStorage.setItem(EXPANDED_KEY, JSON.stringify(next))
  // storage 事件只在其他标签页触发，这里手动通知本标签页订阅者
  expandedCache = null
  expandedListeners.forEach(cb => cb())
}

interface SidebarContentProps {
  chapters: ProjectChapter[]
  pathname: string
  expandedChapters: string[]
  onToggleChapter: (chapterId: string) => void
  onNavigate: () => void
}

function SidebarContent({
  chapters,
  pathname,
  expandedChapters,
  onToggleChapter,
  onNavigate,
}: SidebarContentProps) {
  return (
    <nav className="p-4 space-y-8">
      <div>
        <h3 className="font-medium text-gray-900 dark:text-gray-100 mb-2">
          <Link
            href="/projects"
            className={`inline-flex items-center space-x-3 ${
              pathname === '/projects'
                ? 'text-blue-600 dark:text-blue-400'
                : 'text-gray-600 dark:text-gray-300 hover:text-blue-600 dark:hover:text-blue-400'
            }`}
          >
            <svg
              className="w-5 h-5"
              fill="none"
              viewBox="0 0 24 24"
              stroke="currentColor"
            >
              <path
                strokeLinecap="round"
                strokeLinejoin="round"
                strokeWidth={2}
                d="M3 12l2-2m0 0l7-7 7 7M5 10v10a1 1 0 001 1h3m10-11l2 2m-2-2v10a1 1 0 01-1 1h-3m-6 0a1 1 0 001-1v-4a1 1 0 011-1h2a1 1 0 011 1v4a1 1 0 001 1m-6 0h6"
              />
            </svg>
            <span className="text-base">返回首页</span>
          </Link>
        </h3>
      </div>

      {chapters.map(chapter => (
        <div key={chapter.id}>
          <button
            onClick={() => onToggleChapter(chapter.id)}
            className="w-full flex items-center justify-between font-medium text-gray-900 dark:text-gray-100 mb-2 group hover:text-blue-600 dark:hover:text-blue-400"
          >
            <span>{chapter.title}</span>
            <svg
              className={`w-4 h-4 transition-transform duration-200 ${
                expandedChapters.includes(chapter.id) ? 'transform rotate-180' : ''
              }`}
              fill="none"
              viewBox="0 0 24 24"
              stroke="currentColor"
            >
              <path
                strokeLinecap="round"
                strokeLinejoin="round"
                strokeWidth={2}
                d="M19 9l-7 7-7-7"
              />
            </svg>
          </button>
          {expandedChapters.includes(chapter.id) && (
            <ul className="space-y-2 ml-4">
              {chapter.docs.map(doc => (
                <li key={doc.slug}>
                  <Link
                    href={`/projects/${doc.slug}`}
                    className={`block text-sm ${
                      pathname === `/projects/${doc.slug}`
                        ? 'text-blue-600 dark:text-blue-400'
                        : 'text-gray-600 dark:text-gray-300 hover:text-blue-600 dark:hover:text-blue-400'
                    }`}
                    onClick={onNavigate}
                  >
                    {doc.title}
                  </Link>
                </li>
              ))}
            </ul>
          )}
        </div>
      ))}
    </nav>
  )
}

interface ProjectSidebarProps {
  chapters: ProjectChapter[]
}

export function ProjectSidebar({ chapters }: ProjectSidebarProps) {
  const pathname = usePathname()
  const [isMobileMenuOpen, setIsMobileMenuOpen] = useState(false)
  const expandedChapters = useSyncExternalStore(
    subscribeExpanded,
    getExpandedSnapshot,
    getExpandedServerSnapshot,
  )

  const toggleChapter = (chapterId: string) => {
    const newExpanded = expandedChapters.includes(chapterId)
      ? expandedChapters.filter(id => id !== chapterId)
      : [...expandedChapters, chapterId]
    writeExpanded(newExpanded)
  }

  const sidebarContent = (
    <SidebarContent
      chapters={chapters}
      pathname={pathname}
      expandedChapters={expandedChapters}
      onToggleChapter={toggleChapter}
      onNavigate={() => setIsMobileMenuOpen(false)}
    />
  )

  return (
    <>
      {/* 桌面端侧边栏 */}
      <aside className="hidden md:block w-64 h-[calc(100vh-4rem)] overflow-y-auto sticky top-16 border-r border-gray-200 dark:border-gray-800">
        {sidebarContent}
      </aside>

      {/* 移动端汉堡菜单按钮 */}
      <button
        onClick={() => setIsMobileMenuOpen(!isMobileMenuOpen)}
        className="md:hidden fixed bottom-6 left-6 z-50 p-3 bg-blue-600 text-white rounded-full shadow-lg hover:bg-blue-700 transition-colors"
      >
        {isMobileMenuOpen ? <X size={24} /> : <Menu size={24} />}
      </button>

      {/* 移动端侧边栏 */}
      {isMobileMenuOpen && (
        <div className="md:hidden fixed inset-0 z-40 bg-white dark:bg-gray-900">
          <div className="h-full overflow-y-auto">{sidebarContent}</div>
        </div>
      )}
    </>
  )
}
