local iso = require 'npge.fragment.is_subfragment_of'

return {
    reverse = require 'npge.fragment.reverse',
    has_pos = require 'npge.fragment.has_pos',
    is_subfragment_of = iso, -- FIXME
    subfragment = require 'npge.fragment.subfragment',
    sub = require 'npge.fragment.sub',
    at = require 'npge.fragment.at',
    sequence_to_fragment = require 'npge.fragment.sequence_to_fragment',
    fragment_to_sequence = require 'npge.fragment.fragment_to_sequence',
}