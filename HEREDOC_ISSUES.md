# Here-Document Issues and Alternatives

## Problem
Here-documents (`<<EOF`) consistently fail or get stuck in this environment, causing command execution to hang.

## Examples of Failed Patterns:
```bash
# This pattern fails/hangs:
cat << 'EOF' > file.txt
content here
EOF

# This also fails:
cat > file.txt << 'EOF'
content here
EOF
```

## Working Alternatives:

### 1. Use printf (Recommended)
```bash
printf "content here\nmore content\n" > file.txt
```

### 2. Use echo with newlines
```bash
echo "content here
more content" > file.txt
```

### 3. Use Python scripts
```bash
python3 -c "
with open('file.txt', 'w') as f:
    f.write('content here\nmore content\n')
"
```

### 4. Use separate echo commands
```bash
echo "content here" > file.txt
echo "more content" >> file.txt
```

### 5. Use write_file tool when available
```python
# Use the write_file function directly
```

## Best Practices:
- **Avoid here-docs completely** in this environment
- **Prefer printf** for multi-line content
- **Use Python scripts** for complex content generation
- **Keep commands simple and direct**

## Note:
This appears to be an environment-specific issue with the shell or terminal handling. The alternatives work reliably and should be used instead.
