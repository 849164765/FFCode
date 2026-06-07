#!/usr/bin/env python3
import os
import sys
import pathlib
from pathlib import Path

def convert_to_md(file_path, output_path):
    """将代码文件转换为Markdown格式"""
    try:
        # 读取源文件内容
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # 获取文件扩展名
        ext = file_path.suffix.lstrip('.')
        if ext == '':  # 处理没有扩展名的情况
            ext = 'cpp'
        
        # 构建Markdown内容
        md_content = f"# {file_path.name}\n\n"
        md_content += f"```{ext}\n"
        md_content += content
        md_content += "\n```\n"
        
        # 创建输出目录（如果不存在）
        output_path.parent.mkdir(parents=True, exist_ok=True)
        
        # 写入Markdown文件
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(md_content)
            
        return True
    except Exception as e:
        print(f"转换失败 {file_path}: {str(e)}")
        return False

def main():
    # 定义要处理的文件扩展名
    code_extensions = ['.cpp', '.h', '.hh']
    
    # 获取当前目录
    current_dir = Path.cwd()
    
    # 定义输出目录
    output_base_dir = current_dir / 'src_md'
    
    # 创建输出目录（如果不存在）
    output_base_dir.mkdir(exist_ok=True)
    
    # 遍历当前目录及子目录
    converted_count = 0
    failed_count = 0
    
    for root, _, files in os.walk(current_dir):
        for file in files:
            file_path = Path(root) / file
            
            # 检查文件扩展名
            if file_path.suffix in code_extensions:
                # 计算相对于当前目录的路径
                relative_path = file_path.relative_to(current_dir)
                
                # 构建输出文件路径
                output_file = output_base_dir / relative_path.with_suffix(relative_path.suffix + '.md')
                
                print(f"正在转换: {relative_path}")
                
                # 转换文件
                if convert_to_md(file_path, output_file):
                    converted_count += 1
                else:
                    failed_count += 1
    
    # 打印统计信息
    print("\n转换完成！")
    print(f"成功转换: {converted_count} 个文件")
    print(f"转换失败: {failed_count} 个文件")
    print(f"所有转换后的文件已保存到: {output_base_dir}")

if __name__ == "__main__":
    main()