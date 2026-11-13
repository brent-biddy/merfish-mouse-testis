-- Keep only the figures, for the pptx render only.
--
-- The gfm page is the full document; the deck is its figures. Each figure becomes one slide
-- titled with its `##` heading, so a section that draws several gets several slides with the
-- same title rather than untitled overflow slides. A `#` phase header is kept only when its
-- phase draws at least one figure. Prose, printed output, tables and figure-less sections are
-- dropped. The title slide comes from the metadata and is untouched.
--
-- Registered under the pptx format in the frontmatter, so the gfm render never sees it.
function Pandoc(doc)
  local slides = {}
  local phase = nil     -- the `#` header still waiting for its phase's first figure
  local section = nil   -- the current `##` header

  for _, block in ipairs(doc.blocks) do
    if block.t == "Header" and block.level == 1 then
      phase, section = block, nil
    elseif block.t == "Header" and block.level == 2 then
      section = block
    elseif section then
      block:walk({
        Image = function(image)
          if phase then
            table.insert(slides, phase)
            phase = nil
          end
          image.caption = {}
          table.insert(slides, pandoc.Header(2, section.content))
          table.insert(slides, pandoc.Para({ image }))
        end,
      })
    end
  end

  return pandoc.Pandoc(slides, doc.meta)
end
