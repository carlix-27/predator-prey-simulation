// Ejecutar: node tests/test_explorador.js
const fs = require('fs');
const vm = require('vm');
const assert = require('assert');
const fuente = fs.readFileSync('docs/presentacion_fuente.html', 'utf8');
const inicio = fuente.indexOf('const PARAMS =');
const fin = fuente.indexOf('// =============================================================================\n// Utilidades de dibujo', inicio);
assert(inicio >= 0 && fin > inicio);
vm.runInThisContext(fuente.slice(inicio, fin) + `
  const parametros = { ...PARAMS, K: 11.15 };
  const solucion = integrar(pasoRK2, Y0, 6000, 0.05, parametros);
  if (clasificarRegimen(solucion, 5000, parametros).tipo !== 'transitorio')
    throw new Error('El transitorio lentamente amortiguado debe ser provisional');
  const corta = integrar(pasoRK2, Y0, 1, 0.3, PARAMS);
  if (Math.abs(corta.t[corta.t.length - 1] - 1) > 1e-12) throw new Error('Tiempo final incorrecto');
`);
console.log('Explorador: comprobaciones correctas');
