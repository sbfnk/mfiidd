-- Pandoc/Quarto filter: inline an external file's contents into a CodeBlock
-- whose attributes set `include="path"`. Used to keep session source displays
-- in sync with the actual files in src/ — no copy-paste, no line-number drift.

local function read_file(path)
  local f = io.open(path, "r")
  if not f then
    error("include-code-files: cannot open " .. path)
  end
  local contents = f:read("*all")
  f:close()
  return contents
end

function CodeBlock(block)
  local path = block.attributes["include"]
  if not path then
    return nil
  end
  block.text = read_file(path):gsub("[\r\n]+$", "")
  block.attributes["include"] = nil
  return block
end
