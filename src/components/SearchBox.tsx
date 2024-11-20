import React, { useState } from 'react'

const SearchBox = () => {
  const [searchResults, setSearchResults] = useState([])

  const handleSearch = async (query: string) => {
    const response = await fetch('/api/search')
    const allPosts = await response.json()
    
    const filtered = allPosts.filter(post => 
      post.title.toLowerCase().includes(query.toLowerCase())
    )
    setSearchResults(filtered)
  }

  return (
    <div>
      {/* 添加搜索框组件 */}
    </div>
  )
}

export default SearchBox 