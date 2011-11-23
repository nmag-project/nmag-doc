from sphinx.util.compat import Directive
from sphinx.builders.html import SingleFileHTMLBuilder
from docutils import nodes
from docutils.parsers.rst import directives

class globalindex(nodes.General, nodes.Element):
    pass

def visit_globalindex_node(self, node):
    self.body.append(node['content'])

def depart_globalindex_node(self, node):
    pass

class GlobalIndexDirective(Directive):
    required_arguments = 0
    optional_arguments = 1
    final_argument_whitespace = True
    option_spec = \
      {'maxdepth': directives.nonnegative_int,
       'collapse': directives.flag,
       'titlesonly': directives.flag}

    def run(self):
        node = globalindex('')
        node['maxdepth'] = self.options.get('maxdepth', 2)
        node['collapse'] = 'collapse' in self.options
        node['titlesonly'] = 'titlesonly' in self.options
        return [node]

def process_globalindex_nodes(app, doctree, fromdocname):
    builder = app.builder
    if builder.name != SingleFileHTMLBuilder.name:
        for node in doctree.traverse(globalindex):
            node.parent.remove(node)

    else:
        docname = builder.config.master_doc
        for node in doctree.traverse(globalindex):
            kwargs = dict(maxdepth=node['maxdepth'],
                          collapse=node['collapse'],
                          titles_only=node['titlesonly'])
            rendered_toctree = builder._get_local_toctree(docname, **kwargs)
            node['content'] = rendered_toctree

def setup(app):
    app.add_node(globalindex,
                 html=(visit_globalindex_node, depart_globalindex_node))
    app.add_directive('globalindex', GlobalIndexDirective)
    app.connect('doctree-resolved', process_globalindex_nodes)
