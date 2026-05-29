'use client'

import { useEffect } from 'react'
import { useTheme } from 'next-themes'

/**
 * 自动根据时间设置主题（仅在首次加载且用户未手动切换时生效）
 * 日间 8:00-18:00 使用 light，其余时间使用 dark
 */
export function useAutoTheme() {
  const { setTheme } = useTheme()

  useEffect(() => {
    // next-themes 会在 localStorage 中存储 'theme' key
    // 如果用户已经手动选择过主题，则不自动切换
    if (localStorage.getItem('theme')) return

    const hour = new Date().getHours()
    setTheme(hour >= 8 && hour < 18 ? 'light' : 'dark')
  }, [setTheme])
}
