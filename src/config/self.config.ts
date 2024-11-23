interface SiteConfig {
  // 基础信息
  title: string
  description: string
  author: string
  language: string
  siteUrl: string
  
  // 社交媒体链接
  social: {
    github: string
    email: string
  }
}

const selfConfig: SiteConfig = {
  title: '生物信息文件夹',
  description: '生信工作学习记录',
  author: 'pzweuj',
  language: 'zh-CN',
  siteUrl: 'https://pzweuj.github.io',

  social: {
    github: 'https://github.com/pzweuj',
    email: 'pzweuj@live.com',
  }
}

export default selfConfig
