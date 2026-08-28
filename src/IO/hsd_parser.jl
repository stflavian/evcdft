"""
HSD (Hierarchical Structured Data) Parser for DFTB+ input files.

This module provides a parser for the HSD format used by DFTB+.
The parser supports:
- Key-value pairs: Key = Value
- Blocks: BlockName = { ... }
- Nested blocks
- Arrays: (1, 2, 3) or [1, 2, 3]
- Units: [Angstrom], [Ha], etc.
- Strings in quotes: "value"
- Comments: # comment or // comment
"""

module HSDParser

using ..Constants
using ..Units

# Export the main types and functions
export HSDNode, parse_hsd_file, parse_hsd_string

"""
    Token types for the HSD lexer.
"""
@enum TokenType BEGIN_BLOCK END_BLOCK EQUALS COMMA LPAREN RPAREN LBRACKET RBRACKET STRING NUMBER IDENTIFIER EOF COMMENT

"""
    Represents a token in the input stream.
"""
struct Token
    type::TokenType
    value::String
    line::Int
    column::Int
end

"""
    Represents a node in the HSD tree.
    
    Fields:
    - name: The name/key of this node
    - value: The value (for leaf nodes) or nothing (for blocks)
    - children: Child nodes (for blocks)
    - units: Units specification (e.g., "Angstrom", "Ha")
    - is_array: Whether this node represents an array
"""
mutable struct HSDNode
    name::String
    value::Union{Nothing, String, Float64, Int, Vector{Float64}, Vector{Int}, Vector{String}}
    children::Dict{String, HSDNode}
    units::Union{Nothing, String}
    is_array::Bool
    
    function HSDNode(name::String)
        new(name, nothing, Dict{String, HSDNode}(), nothing, false)
    end
    
    function HSDNode(name::String, value)
        new(name, value, Dict{String, HSDNode}(), nothing, false)
    end
end

"""
    Lexer: Convert input string to tokens.
    
    Args:
    - input: The input string to tokenize
    
    Returns:
    - Vector of tokens
"""
function tokenize(input::String)::Vector{Token}
    tokens = Token[]
    pos = 1
    line = 1
    column = 1
    len = length(input)
    
    while pos <= len
        char = input[pos]
        
        # Skip whitespace
        if isspace(char)
            if char == '\n'
                line += 1
                column = 1
            else
                column += 1
            end
            pos += 1
            continue
        end
        
        # Handle comments
        if char == '#' || (char == '/' && pos < len && input[pos+1] == '/')
            # Skip to end of line
            while pos <= len && input[pos] != '\n'
                pos += 1
            end
            if pos <= len && input[pos] == '\n'
                line += 1
                column = 1
                pos += 1
            end
            continue
        end
        
        # Handle identifiers and keywords
        if isletter(char) || char == '_'
            start = pos
            while pos <= len && (isletter(input[pos]) || isdigit(input[pos]) || input[pos] == '_' || input[pos] == '-')
                pos += 1
                column += 1
            end
            value = input[start:pos-1]
            
            # Check for special keywords
            if value == "{" 
                push!(tokens, Token(BEGIN_BLOCK, "{", line, column - length(value)))
            elseif value == "}"
                push!(tokens, Token(END_BLOCK, "}", line, column - length(value)))
            elseif value == "=" 
                push!(tokens, Token(EQUALS, "=", line, column - length(value)))
            else
                push!(tokens, Token(IDENTIFIER, value, line, column - length(value)))
            end
            continue
        end
        
        # Handle numbers
        if isdigit(char) || (char == '-' && pos < len && isdigit(input[pos+1]))
            start = pos
            has_decimal = false
            has_exponent = false
            
            if char == '-'
                pos += 1
                column += 1
            end
            
            while pos <= len
                c = input[pos]
                if isdigit(c)
                    pos += 1
                    column += 1
                elseif c == '.' && !has_decimal
                    has_decimal = true
                    pos += 1
                    column += 1
                elseif (c == 'e' || c == 'E') && !has_exponent
                    has_exponent = true
                    pos += 1
                    column += 1
                    if pos <= len && (input[pos] == '+' || input[pos] == '-')
                        pos += 1
                        column += 1
                    end
                else
                    break
                end
            end
            
            value = input[start:pos-1]
            push!(tokens, Token(NUMBER, value, line, column - length(value)))
            continue
        end
        
        # Handle single-character tokens
        if char == '='
            push!(tokens, Token(EQUALS, "=", line, column))
            pos += 1
            column += 1
        elseif char == ','
            push!(tokens, Token(COMMA, ",", line, column))
            pos += 1
            column += 1
        elseif char == '('
            push!(tokens, Token(LPAREN, "(", line, column))
            pos += 1
            column += 1
        elseif char == ')'
            push!(tokens, Token(RPAREN, ")", line, column))
            pos += 1
            column += 1
        elseif char == '['
            push!(tokens, Token(LBRACKET, "[", line, column))
            pos += 1
            column += 1
        elseif char == ']'
            push!(tokens, Token(RBRACKET, "]", line, column))
            pos += 1
            column += 1
        elseif char == '{'
            push!(tokens, Token(BEGIN_BLOCK, "{", line, column))
            pos += 1
            column += 1
        elseif char == '}'
            push!(tokens, Token(END_BLOCK, "}", line, column))
            pos += 1
            column += 1
        elseif char == '"'
            # Handle quoted string
            start = pos
            pos += 1
            column += 1
            while pos <= len && input[pos] != '"'
                if input[pos] == '\\' && pos < len
                    pos += 2
                    column += 2
                else
                    pos += 1
                    column += 1
                end
            end
            if pos <= len && input[pos] == '"'
                pos += 1
                column += 1
            end
            value = input[start+1:pos-2]  # Remove quotes
            push!(tokens, Token(STRING, value, line, column - length(value) - 2))
            continue
        else
            # Unknown character, skip
            pos += 1
            column += 1
        end
    end
    
    push!(tokens, Token(EOF, "", line, column))
    return tokens
end

"""
    Parser: Convert tokens to HSDNode tree.
    
    Uses recursive descent parsing.
    
    Args:
    - tokens: Vector of tokens from lexer
    - pos: Current position in token stream (1-indexed)
    
    Returns:
    - root: Root HSDNode
    - new_pos: Updated position
"""
function parse_tokens(tokens::Vector{Token}, pos::Int=1)::Tuple{HSDNode, Int}
    root = HSDNode("root")
    return parse_block(tokens, pos, root)
end

"""
    Parse a block (enclosed in { } or at top level).
"""
function parse_block(tokens::Vector{Token}, pos::Int, parent::HSDNode)::Tuple{HSDNode, Int}
    while pos <= length(tokens)
        token = tokens[pos]
        
        if token.type == END_BLOCK || token.type == EOF
            break
        end
        
        if token.type == COMMENT
            pos += 1
            continue
        end
        
        # Parse a key-value pair or block
        if token.type == IDENTIFIER
            key = token.value
            pos += 1
            
            # Check if next token is EQUALS
            if pos <= length(tokens) && tokens[pos].type == EQUALS
                pos += 1  # Skip EQUALS
                
                # Parse the value
                value, pos = parse_value(tokens, pos)
                
                # Check for units
                units = nothing
                if pos <= length(tokens) && tokens[pos].type == LBRACKET
                    pos += 1  # Skip [
                    if pos <= length(tokens) && tokens[pos].type == IDENTIFIER
                        units = tokens[pos].value
                        pos += 1
                    end
                    if pos <= length(tokens) && tokens[pos].type == RBRACKET
                        pos += 1
                    end
                end
                
                # Create node with value
                node = HSDNode(key, value)
                node.units = units
                parent.children[key] = node
            elseif pos <= length(tokens) && tokens[pos].type == BEGIN_BLOCK
                pos += 1  # Skip {
                
                # Create block node
                node = HSDNode(key)
                parent.children[key] = node
                
                # Parse the block contents
                node, pos = parse_block(tokens, pos, node)
                
                # Skip closing }
                if pos <= length(tokens) && tokens[pos].type == END_BLOCK
                    pos += 1
                end
            else
                # Just an identifier without value - treat as flag
                node = HSDNode(key, true)
                parent.children[key] = node
            end
        elseif token.type == BEGIN_BLOCK
            # Anonymous block - should not happen in valid HSD
            pos += 1
        elseif token.type == NUMBER
            # Handle bare numbers in blocks (e.g., matrix rows without keys)
            # Collect all consecutive numbers on this line as a vector
            row_values = Float64[]
            while pos <= length(tokens) && (tokens[pos].type == NUMBER || tokens[pos].type == COMMA)
                if tokens[pos].type == NUMBER
                    push!(row_values, parse(Float64, tokens[pos].value))
                end
                pos += 1
            end
            # Store as an anonymous child with index-based key
            key = string(length(parent.children) + 1)
            node = HSDNode(key, row_values)
            parent.children[key] = node
        else
            pos += 1
        end
    end
    
    return parent, pos
end

"""
    Parse a value (number, string, array, or identifier).
"""
function parse_value(tokens::Vector{Token}, pos::Int)::Tuple{Union{Nothing, String, Float64, Int, Vector}, Int}
    token = tokens[pos]
    
    if token.type == NUMBER
        # Try to parse as Int first, then Float64
        try
            value = parse(Int, token.value)
            return value, pos + 1
        catch
            try
                value = parse(Float64, token.value)
                return value, pos + 1
            catch
                return token.value, pos + 1
            end
        end
    elseif token.type == STRING
        return token.value, pos + 1
    elseif token.type == IDENTIFIER
        # Could be a boolean or string
        if lowercase(token.value) == "yes" || lowercase(token.value) == "true"
            return true, pos + 1
        elseif lowercase(token.value) == "no" || lowercase(token.value) == "false"
            return false, pos + 1
        else
            return token.value, pos + 1
        end
    elseif token.type == LPAREN
        # Parse array
        pos += 1  # Skip (
        elements = []
        while pos <= length(tokens) && tokens[pos].type != RPAREN
            if tokens[pos].type == COMMA
                pos += 1
                continue
            end
            value, pos = parse_value(tokens, pos)
            push!(elements, value)
        end
        if pos <= length(tokens) && tokens[pos].type == RPAREN
            pos += 1
        end
        return elements, pos
    elseif token.type == LBRACKET
        # This might be a unit specification, handled elsewhere
        return nothing, pos
    elseif token.type == BEGIN_BLOCK
        # Parse { ... } as array literal
        pos += 1  # Skip {
        elements = []
        while pos <= length(tokens) && tokens[pos].type != END_BLOCK
            if tokens[pos].type == COMMA
                pos += 1
                continue
            end
            value, pos = parse_value(tokens, pos)
            if value !== nothing
                push!(elements, value)
            end
        end
        if pos <= length(tokens) && tokens[pos].type == END_BLOCK
            pos += 1  # Skip }
        end
        return elements, pos
    else
        return nothing, pos + 1
    end
end

"""
    Parse an HSD file and return the root node.
    
    Args:
    - filename: Path to the HSD file
    
    Returns:
    - Root HSDNode
"""
function parse_hsd_file(filename::String)::HSDNode
    content = read(filename, String)
    return parse_hsd_string(content)
end

"""
    Parse an HSD string and return the root node.
    
    Args:
    - content: HSD content as a string
    
    Returns:
    - Root HSDNode
"""
function parse_hsd_string(content::String)::HSDNode
    tokens = tokenize(content)
    root, _ = parse_tokens(tokens)
    return root
end

"""
    Get a child node by path (e.g., "Geometry.LatticeVectors").
    
    Args:
    - node: Parent node
    - path: Dot-separated path to child
    
    Returns:
    - Child node or nothing if not found
"""
function get_node(node::HSDNode, path::String)::Union{HSDNode, Nothing}
    parts = split(path, '.')
    current = node
    for part in parts
        if haskey(current.children, part)
            current = current.children[part]
        else
            return nothing
        end
    end
    return current
end

"""
    Get the value of a node, converting units if specified.
    
    Args:
    - node: HSDNode to get value from
    
    Returns:
    - Value with units converted (if applicable)
"""
function get_value(node::HSDNode)
    value = node.value
    
    if value === nothing
        return nothing
    end
    
    if node.units === nothing
        return value
    end
    
    # Handle unit conversion
    if node.units == "Angstrom" || node.units == "angstrom"
        if value isa Real
            return value * angstrom_to_bohr
        elseif value isa Vector
            return value .* angstrom_to_bohr
        end
    elseif node.units == "Bohr" || node.units == "bohr"
        # Already in Bohr, no conversion needed
        return value
    elseif node.units == "Ha" || node.units == "hartree"
        # Already in Hartree
        return value
    elseif node.units == "eV"
        if value isa Real
            return value * ev_to_hartree
        elseif value isa Vector
            return value .* ev_to_hartree
        end
    elseif node.units == "nm" || node.units == "nanometer"
        if value isa Real
            return value * 10.0 * angstrom_to_bohr  # 1 nm = 10 Å
        elseif value isa Vector
            return value .* 10.0 .* angstrom_to_bohr
        end
    end
    
    return value
end

"""
    Pretty print an HSDNode tree (for debugging).
"""
function print_node(node::HSDNode, indent::Int=0)
    prefix = "  "^indent
    
    if node.value !== nothing
        units_str = node.units === nothing ? "" : " [$node.units]"
        if node.is_array
            println(prefix * "$node.name = $node.value$units_str")
        else
            println(prefix * "$node.name = $node.value$units_str")
        end
    else
        println(prefix * "$node.name = {")
        for (key, child) in node.children
            print_node(child, indent + 1)
        end
        println(prefix * "}")
    end
end

end # module HSDParser
