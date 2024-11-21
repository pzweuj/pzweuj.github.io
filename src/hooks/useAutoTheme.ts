'use client'

import { useEffect, useState } from 'react'
import { useTheme } from 'next-themes'

export function useAutoTheme() {
  const { theme, setTheme } = useTheme()
  const [mounted, setMounted] = useState(false)

  useEffect(() => {
    setMounted(true)
  }, [])

  useEffect(() => {
    if (typeof window === 'undefined') return

    // 获取保存的用户偏好
    const savedTheme = localStorage.getItem('user-theme-preference')
    
    if (savedTheme) {
      setTheme(savedTheme)
      return
    }

    // 检查当前时间
    const hour = new Date().getHours()
    const shouldBeDark = hour < 8 || hour >= 18
    
    setTheme(shouldBeDark ? 'dark' : 'light')
  }, [setTheme])

  const toggleTheme = () => {
    if (typeof window === 'undefined') return

    const newTheme = theme === 'dark' ? 'light' : 'dark'
    setTheme(newTheme)
    localStorage.setItem('user-theme-preference', newTheme)
  }

  // 如果组件未挂载，返回null或初始值
  if (!mounted) {
    return { theme: undefined, toggleTheme }
  }

  return { theme, toggleTheme }
} 