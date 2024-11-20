'use client'

import { useEffect } from 'react'
import { useTheme } from 'next-themes'

'use client'

import { useEffect } from 'react'
import { useTheme } from 'next-themes'

export function useAutoTheme() {
  const { theme, setTheme } = useTheme()

  useEffect(() => {
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
    const newTheme = theme === 'dark' ? 'light' : 'dark'
    setTheme(newTheme)
    localStorage.setItem('user-theme-preference', newTheme)
  }

  return { theme, toggleTheme }
} 