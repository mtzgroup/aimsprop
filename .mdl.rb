# Enable all rules by default
all

# Allow any line length
exclude_rule 'MD013'

# Allow ordered lists
rule 'MD029', :style => 'ordered'

# Allow "?" in headers
rule 'MD026', :punctuation => '.,;:!'

# Allow multiple headers to be the same, as required for changelogs
# https://github.com/markdownlint/markdownlint/blob/master/docs/RULES.md#md024---multiple-headers-with-the-same-content
rule 'MD024', :allow_different_nesting => true