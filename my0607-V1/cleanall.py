import os
import re

def remove_comments_before_marker():
    # 遍历当前文件夹及其子文件夹
    for root, dirs, files in os.walk('.'):
        for file in files:
            # 只处理代码文件
            if file.lower().endswith(('.cpp', '.cc', '.cxx', '.c', '.h', '.hpp', '.hh', '.hxx', '.py', '.js', '.ts', '.java', '.cs', '.php', '.rb', '.go', '.rs', '.swift', '.kt', '.scala')):
                file_path = os.path.join(root, file)
                
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        content = f.read()
                    
                    # 找到文件名在内容中的首次出现位置
                    file_name_pattern = r'//\s*' + re.escape(file) + r'\s*$'
                    matches = list(re.finditer(file_name_pattern, content, re.MULTILINE))
                    
                    if matches:
                        # 获取首次出现的位置
                        first_occurrence_pos = matches[0].start()
                        
                        # 提取标记前的内容
                        content_before_marker = content[:first_occurrence_pos]
                        content_after_marker = content[first_occurrence_pos:]
                        
                        # 删除标记前的注释和空行
                        # 删除单行注释 (// ...)
                        cleaned_content_before = re.sub(r'^\s*//.*$', '', content_before_marker, flags=re.MULTILINE)
                        
                        # 删除多行注释 (/* ... */)
                        cleaned_content_before = re.sub(r'/\*.*?\*/', '', cleaned_content_before, flags=re.DOTALL)
                        
                        # 删除空行和只包含空白字符的行
                        lines = cleaned_content_before.splitlines()
                        non_empty_lines = [line for line in lines if line.strip() != '']
                        cleaned_content_before = '\n'.join(non_empty_lines)
                        
                        # 重新组合内容
                        new_content = cleaned_content_before + content_after_marker
                        
                        # 再次清理可能产生的多余空行
                        new_content = re.sub(r'\n\s*\n', '\n\n', new_content)  # 将多个空行合并为两个
                        new_content = new_content.strip()  # 去除首尾空白
                        
                        # 写回文件
                        with open(file_path, 'w', encoding='utf-8') as f:
                            f.write(new_content)
                        
                        print(f"已处理文件: {file_path}")
                    else:
                        print(f"文件 {file_path} 中未找到标记 '// {file}'")
                        
                except UnicodeDecodeError:
                    # 如果不是文本文件，跳过
                    print(f"跳过非文本文件: {file_path}")
                except Exception as e:
                    print(f"处理文件 {file_path} 时出错: {str(e)}")

def remove_comments_before_custom_marker(marker):
    """
    通用函数：删除指定标记前的所有注释和空行
    :param marker: 要查找的标记，例如 "// basic_geometry.hh"
    """
    # 遍历当前文件夹及其子文件夹
    for root, dirs, files in os.walk('.'):
        for file in files:
            # 只处理代码文件
            if file.lower().endswith(('.cpp', '.cc', '.cxx', '.c', '.h', '.hpp', '.hh', '.hxx', '.py', '.js', '.ts', '.java', '.cs', '.php', '.rb', '.go', '.rs', '.swift', '.kt', '.scala')):
                file_path = os.path.join(root, file)
                
                try:
                    with open(file_path, 'r', encoding='utf-8') as f:
                        content = f.read()
                    
                    # 查找标记首次出现的位置
                    marker_pos = content.find(marker)
                    
                    if marker_pos != -1:
                        # 提取标记前的内容
                        content_before_marker = content[:marker_pos]
                        content_after_marker = content[marker_pos:]
                        
                        # 删除标记前的注释和空行
                        # 删除单行注释 (// ...)
                        cleaned_content_before = re.sub(r'^\s*//.*$', '', content_before_marker, flags=re.MULTILINE)
                        
                        # 删除多行注释 (/* ... */)
                        cleaned_content_before = re.sub(r'/\*.*?\*/', '', cleaned_content_before, flags=re.DOTALL)
                        
                        # 删除空行和只包含空白字符的行
                        lines = cleaned_content_before.splitlines()
                        non_empty_lines = [line for line in lines if line.strip() != '']
                        cleaned_content_before = '\n'.join(non_empty_lines)
                        
                        # 重新组合内容
                        new_content = cleaned_content_before + content_after_marker
                        
                        # 再次清理可能产生的多余空行
                        new_content = re.sub(r'\n\s*\n', '\n\n', new_content)  # 将多个空行合并为两个
                        new_content = new_content.strip()  # 去除首尾空白
                        
                        # 写回文件
                        with open(file_path, 'w', encoding='utf-8') as f:
                            f.write(new_content)
                        
                        print(f"已处理文件: {file_path}，标记: {marker}")
                    else:
                        print(f"文件 {file_path} 中未找到标记: {marker}")
                        
                except UnicodeDecodeError:
                    # 如果不是文本文件，跳过
                    print(f"跳过非文本文件: {file_path}")
                except Exception as e:
                    print(f"处理文件 {file_path} 时出错: {str(e)}")

def interactive_mode():
    """交互模式，让用户选择操作方式"""
    print("请选择操作模式:")
    print("1. 删除文件名标记前的注释 (例如 '// filename.ext')")
    print("2. 删除自定义标记前的注释")
    
    choice = input("请输入选择 (1 或 2): ").strip()
    
    if choice == '1':
        print("开始删除文件名标记前的注释...")
        remove_comments_before_marker()
    elif choice == '2':
        marker = input("请输入要查找的标记 (例如 '// basic_geometry.hh'): ").strip()
        print(f"开始删除标记 '{marker}' 前的注释...")
        remove_comments_before_custom_marker(marker)
    else:
        print("无效选择，退出程序。")

if __name__ == "__main__":
    print("开始处理文件...")
    interactive_mode()
    print("处理完成！")



