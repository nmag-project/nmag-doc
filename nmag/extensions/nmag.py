from sphinx.util.compat import Directive
from sphinx.builders.html import SingleFileHTMLBuilder
from docutils import nodes
from docutils.parsers.rst import directives

class globalindex(nodes.General, nodes.Element):
    pass

def visit_globalindex_node(self, node):
    self.body.append(node['my_content'])

def depart_globalindex_node(self, node):
    pass

class GlobalIndexDirective(Directive):
    has_content = True
    required_arguments = 0
    optional_arguments = 0
    final_argument_whitespace = True
    option_spec = \
      {'maxdepth': directives.nonnegative_int,
       'collapse': directives.flag}

    def run(self):
        node = globalindex('')
        node['maxdepth'] = self.options.get('maxdepth', 2)
        node['collapse'] = 'collapse' in self.options
        print node['my_content']
        raw_input()
        return [node]

def process_globalindex_nodes(app, doctree, fromdocname):
    builder = app.builder
    if builder.name != 'nmaghtml':
        for node in doctree.traverse(globalindex):
            node.parent.remove(node)

    else:
        docname = builder.config.master_doc
        for node in doctree.traverse(globalindex):
            kwargs = dict(maxdepth=2, collapse=False)
            rendered_toctree = builder._get_local_toctree(docname, **kwargs)
            node['my_content'] = rendered_toctree


class NmagHTMLBuilder(SingleFileHTMLBuilder):
    name = 'nmaghtml'

    def finish(self):
        self.write_genindex()
        SingleFileHTMLBuilder.finish(self)

    def get_doc_context(self, *args):
        ctx = SingleFileHTMLBuilder.get_doc_context(self, *args)
        return ctx

    def write(self, *args):

        x = SingleFileHTMLBuilder.write(self, *args)

        return x

    def _get_local_toctree(self, *args, **kwargs):
         x = SingleFileHTMLBuilder._get_local_toctree(self, *args, **kwargs)
         return x

    def write_genindex(self):
        genindex = self.env.create_index(self)
        indexcounts = []
        for _, entries in genindex:
            indexcounts.append(sum(1 + len(subitems)
                                   for _, (_, subitems) in entries))
        genindexcontext = \
          dict(genindexentries=genindex,
               genindexcounts=indexcounts,
               split_index=False)
        self.handle_page('genindex', genindexcontext, 'genindex.html')

def setup(app):
    app.add_builder(NmagHTMLBuilder)

    app.add_node(globalindex,
                 html=(visit_globalindex_node, depart_globalindex_node),
                 latex=(visit_globalindex_node, depart_globalindex_node),
                 text=(visit_globalindex_node, depart_globalindex_node))

    app.add_directive('globalindex', GlobalIndexDirective)
    app.connect('doctree-resolved', process_globalindex_nodes)
    #app.connect('env-purge-doc', purge_todos)
