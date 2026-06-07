import os
import re

def remove_copyright_content():
    # 定义要删除的版权信息内容
    copyright_content = """"""
    
    # 遍历当前文件夹及其子文件夹
    for root, dirs, files in os.walk('.'):
        for file in files:
            file_path = os.path.join(root, file)
            
            # 检查文件是否为文本文件
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.read()
                
                # 检查文件内容是否包含版权信息
                if copyright_content in content:
                    # 删除版权信息
                    new_content = content.replace(copyright_content, '').strip()
                    
                    # 如果删除后内容为空，可以选择删除整个文件或保留为空
                    with open(file_path, 'w', encoding='utf-8') as f:
                        f.write(new_content)
                    
                    print(f"已从文件 {file_path} 中删除版权信息")
                    
            except UnicodeDecodeError:
                # 如果不是文本文件，跳过
                print(f"跳过非文本文件: {file_path}")
            except Exception as e:
                print(f"处理文件 {file_path} 时出错: {str(e)}")

if __name__ == "__main__":
    print("开始删除包含指定版权信息的文件内容...")
    remove_copyright_content()
    print("处理完成！")