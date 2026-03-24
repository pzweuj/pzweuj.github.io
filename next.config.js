/** @type {import('next').NextConfig} */
const nextConfig = {
  output: 'export',
  images: {
    unoptimized: true,
  },
  // 并行构建优化
  experimental: {
    // 使用 worker 线程并行生成静态页面
    workerThreads: true,
    // 限制并发数，避免内存溢出 (GitHub Actions 默认 2 核)
    cpus: 2,
  },
}

module.exports = nextConfig 