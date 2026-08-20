# Introdução ao método dos elementos de contorno

Notas de aula sobre o Método dos Elementos de Contorno (BEM), em [Typst](https://typst.app/).  
Curso de **30 h**. Código de referência: [`BEM_gmsh`](https://github.com/l-s-campos/BEM_gmsh).

**Site (multipágina):** https://l-s-campos.github.io/Introdu-o-ao-m-todo-dos-elementos-de-contorno/

## Capítulos (trilha principal)

1. [Apresentação](docs/apresentacao.html)
2. [Glossário e notação](docs/glossario.html)
3. [Interpolação](docs/interpolacao.html)
4. [Equações diferenciais](docs/equacoes-diferenciais.html)
5. [Indo para 2D](docs/indo-para-2d.html)
6. [Laplace 2D](docs/laplace-2d.html)
7. [Poisson 2D](docs/poisson-2d.html)
8. [Elasticidade 2D](docs/elasticidade-2d.html)
9. [Indo para 3D](docs/indo-para-3d.html)
10. [Trabalhos finais](docs/trabalhos-finais.html)

## Extra e apêndice

- [Viga de Euler](docs/viga-euler.html) (aprofundamento opcional)
- [Contato half-space](docs/contato-halfspace.html) (opcional)
- [Extra / medidas de erro (apêndice)](docs/extra.html)

## Compilar

```bash
# PDF completo
typst compile main.typ

# Site multipágina (HTML)
node build-site.mjs
```

## Conversão

Markdown → Typst com [`markdown2typst`](https://github.com/Mapaor/markdown2typst) + [`tex2typst`](https://github.com/qwinsi/tex2typst).

## Links

- [Apostila](https://1drv.ms/b/s!AmfyGvdmTYong45aJ5g2TBxKCkygcQ?e=DJQ9oC)
- [Entregas](https://forms.gle/7gKy3k1TqHaCUkVD8)
- [Playlist YouTube](https://www.youtube.com/playlist?list=PLajnQa6HBzEIrJXXrQfUAeYwdMk1ygAhr)
- [Código BEM_gmsh](https://github.com/l-s-campos/BEM_gmsh)
