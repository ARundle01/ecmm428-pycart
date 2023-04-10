<head>
    <script>
    MathJax = {
      tex: {
        inlineMath: [['\\\\(', '\\\\)'], ['$', '$']]
      }
    };
    </script>
    <script id="MathJax-script" async
      src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js">
    </script>
</head> 

::: pycart.cartogram.Cartogram
    handler: python
    options:
        members:
            - __init__
            - non_contiguous
            - dorling
        show_root_heading: true
        show_source: false
---

## Helpers
### `_paired_distance(X, Y)` 
::: pycart.cartogram._paired_distance
    handler: python
    options:
        show_source: false

---

### `_repel(neighbour, focal, xrepel, yrepel)`
::: pycart.cartogram._repel
    handler: python
    options:
        show_source: false

---

### `_attract(nb, borders, idx, focal, perimeter, xattract, yattract)`
::: pycart.cartogram._attract
    handler: python
    options:
        members:
            - _attract
        show_source: false

---

# References
<a href="https://www.tandfonline.com/doi/epdf/10.1111/j.0033-0124.1976.00371.x?needAccess=true" target="_blank">[1]</a>
J. M. Olson, "Noncontiguous Area Cartograms," in *The Professional Geographer*, vol. 28, no. 4, 
pp. 371-380, 1976.

<a href="https://www.dannydorling.org/wp-content/files/dannydorling_publication_id1448.pdf" target="_blank">[2]</a>
D. Dorling, "Area Cartograms: Their Use and Creation," in *Concepts and Techniques in Modern Geography*, no. 59
Environmental Publications, University of East Anglia, pp. 32-36, 1996.

---