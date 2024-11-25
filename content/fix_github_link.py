import os
import re

def fix_github_links(folder_path):
    """
    遍历指定文件夹中的所有md文件，替换githubusercontent链接
    
    Args:
        folder_path: 要搜索的文件夹路径
    """
    # githubusercontent链接的正则表达式模式
    old_pattern = 'https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/master'
    new_pattern = 'https://raw.githubusercontent.com/pzweuj/pzweuj.github.io/refs/heads/master'
    
    # 遍历文件夹
    for root, dirs, files in os.walk(folder_path):
        for file in files:
            if file.endswith('.md'):
                file_path = os.path.join(root, file)
                print(f"处理文件: {file_path}")
                try:
                    # 读取文件内容
                    with open(file_path, 'r', encoding='utf-8') as f:
                        content = f.read()
                    
                    # 查找并替换链接
                    new_content = re.sub(old_pattern, new_pattern, content)
                    
                    # 如果内容有变化，写回文件
                    if new_content != content:
                        with open(file_path, 'w', encoding='utf-8') as f:
                            f.write(new_content)
                        print(f"已更新文件: {file_path}")
                
                except Exception as e:
                    print(f"处理文件 {file_path} 时出错: {str(e)}")

# 使用示例
fix_github_links("posts")
