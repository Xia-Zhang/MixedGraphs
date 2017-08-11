HTMLWidgets.widget({

  name: 'Cytoscapejs',

  type: 'output',

  factory: function(el, width, height) {

    return {

      renderValue: function(x) {
        for (var node of x.nodes) {
          node.data.color = node.data.color.substring(0, 7);
        }
        
        cytoscape({
          container: el,

          layout: {
            name: 'cose',
            padding: 10
          },

          style: [
            {
              selector: 'node',
              style: {
                'background-color': 'data(color)',
                'content': 'data(name)',
                'text-valign': 'center'
              }
            },
            {
              selector: 'edge',
              style: {
                'curve-style': 'bezier',
                'target-arrow-shape': function(ele) { return ele.data('directed') ? 'triangle' : 'none' },
                'source-arrow-shape': 'none',
                'width': function(ele) { return ele.data('weight') ? Math.abs(ele.data('weight') * 10) : 5 },
                'line-color': '#aaaaaa',
                'shape': 'roundrectangle',
                'background-color': '#ccc',
                'background-opacity': 0.333,
              }
            }
          ],

          elements: x,

          ready: function(){
            window.cy = this;
          }
        });
      },

      resize: function(width, height) {
        
        if (instance.cytoscape)
            instance.cytoscape.resize();        

      }   

    };
  }
});