HTMLWidgets.widget({

  name: 'Cytoscapejs',

  type: 'output',

  factory: function(el, width, height) {

    return {

      renderValue: function(x) {
        for (var node of x.nodes) {
          node.data.color = node.data.color.substring(0, 7);
        }
        for (var edge of x.edges) {
          if (edge.data.weight) {
            edge.data.renderWeight = Math.abs(edge.data.weight * 10);  
          } else {
            edge.data.renderWeight = 5;
          }
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
                'target-arrow-shape': 'triangle',
                'source-arrow-shape': 'circle',                
                'width': 'data(renderWeight)',
                'line-color': '#aaaaaa'
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